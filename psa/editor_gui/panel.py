"""Shared panel of the editor: metadata, gate, robustness, 3D preview, fidelity vs. offset, save / send."""
from __future__ import annotations

import os
import traceback

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QTimer, Signal
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
                               QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget)

from psa.design import metrics
from psa.design.conversion import ps_to_curve
from psa.design.pulse import PulseSequence
from psa.editor_gui.widgets import BAD, GOOD, MUTED_FG, Mpl3D, make_plot

ROBUST_THRESHOLD = 1e-6       # normalized |A_k| / T^k below this counts as cancelled


def fmt_vec(v, digits=3):
    return "(" + ", ".join(f"{x:+.{digits}f}" for x in v) + ")"


class SummaryPanel(QWidget):
    metadataChanged = Signal(float, float)       # dt [s], umax [Hz]
    sendRequested = Signal(str)                  # saved file path

    def __init__(self, dt=0.5e-6, umax=1e4, save_path="", parent=None):
        super().__init__(parent)
        self.pulse = None
        self.setMinimumWidth(380)
        layout = QVBoxLayout(self)

        meta = QGroupBox("Metadata")
        form = QFormLayout(meta)
        self.dt_spin = QDoubleSpinBox(); self.dt_spin.setRange(0.001, 1e6); self.dt_spin.setDecimals(4)
        self.dt_spin.setSuffix(" µs"); self.dt_spin.setValue(dt * 1e6)
        self.umax_spin = QDoubleSpinBox(); self.umax_spin.setRange(0.001, 1e6); self.umax_spin.setDecimals(4)
        self.umax_spin.setSuffix(" kHz"); self.umax_spin.setValue(umax / 1e3)
        form.addRow("sub-pulse duration:", self.dt_spin)
        form.addRow("max. amplitude:", self.umax_spin)
        self.dt_spin.editingFinished.connect(self._metadata_edited)
        self.umax_spin.editingFinished.connect(self._metadata_edited)
        layout.addWidget(meta)

        res = QGroupBox("Resulting pulse")
        rform = QFormLayout(res)
        self.labels = {}
        for key, text in (("info", "pulse:"), ("amp", "max. amplitude:"), ("gate", "rotation angle:"),
                          ("axis", "rotation axis:"), ("A1", "1st order |A1|/T:"), ("A2", "2nd order |A2|/T²:"),
                          ("A3", "3rd order |A3|/T³:"), ("B1", "B1 scaling (1st):")):
            lab = QLabel("–")
            self.labels[key] = lab
            rform.addRow(text, lab)
        layout.addWidget(res)

        self.view3d = Mpl3D()
        self.view3d.setMinimumHeight(240)
        layout.addWidget(self.view3d, stretch=2)

        fid = QGroupBox("Fidelity vs. offset")
        fid.setToolTip("gate fidelity (|tr U(0)^† U(offset)|² + 2)/6 relative to the on-resonance gate")
        fl = QVBoxLayout(fid)
        row = QHBoxLayout()
        self.fid_check = QCheckBox("live"); self.fid_check.setChecked(True)
        self.fid_range = QDoubleSpinBox(); self.fid_range.setRange(0.1, 1000); self.fid_range.setValue(20)
        self.fid_range.setSuffix(" kHz"); self.fid_range.setToolTip("total offset range (centred on 0)")
        self.fid_log = QCheckBox("log10(1-F)"); self.fid_log.setChecked(True)
        for w in (self.fid_check, self.fid_range, self.fid_log):
            row.addWidget(w)
        fl.addLayout(row)
        self.fid_plot = make_plot("log10(1-F)", "offset [kHz]")
        self.fid_plot.setMinimumHeight(150)
        self.fid_curve = self.fid_plot.plot(pen=pg.mkPen("#38bdf8", width=2))
        fl.addWidget(self.fid_plot)
        for w in (self.fid_check, self.fid_log):
            w.toggled.connect(self.schedule)
        self.fid_range.editingFinished.connect(self.schedule)
        layout.addWidget(fid, stretch=1)

        out = QGroupBox("Output")
        ol = QVBoxLayout(out)
        prow = QHBoxLayout()
        self.path_edit = QLineEdit(save_path)
        self.path_edit.setToolTip("File to save; the extension (.csv or .bruker) decides the format")
        browse = QPushButton("…"); browse.setFixedWidth(30); browse.clicked.connect(self._browse)
        prow.addWidget(self.path_edit); prow.addWidget(browse)
        ol.addLayout(prow)
        frow = QHBoxLayout()
        frow.addWidget(QLabel("format:"))
        self.format_box = QComboBox()
        self.format_box.addItems(["CSV (.csv)", "Bruker shape (.bruker)"])
        self.format_box.setToolTip("CSV: lines 'amplitude %, phase °'.  Bruker: JCAMP shape file with "
                                   "timestep and nominal RF in the header, phases wrapped to [0, 360).")
        frow.addWidget(self.format_box, stretch=1)
        ol.addLayout(frow)
        self._sync_format_from_path()
        self.format_box.currentIndexChanged.connect(self._format_changed)
        self.path_edit.textEdited.connect(self._sync_format_from_path)
        brow = QHBoxLayout()
        self.save_button = QPushButton("Save")
        self.send_button = QPushButton("Save + plot in main window")
        self.save_button.clicked.connect(self.save)
        self.send_button.clicked.connect(self.save_and_send)
        brow.addWidget(self.save_button); brow.addWidget(self.send_button)
        ol.addLayout(brow)
        self.status = QLabel(""); self.status.setWordWrap(True)
        self.status.setStyleSheet(f"color: {MUTED_FG};")
        ol.addWidget(self.status)
        layout.addWidget(out)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(200)
        self._timer.timeout.connect(self.recompute)

    # ------------------------------------------------------------ metadata ---
    @property
    def dt(self):
        return self.dt_spin.value() * 1e-6

    @property
    def umax(self):
        return self.umax_spin.value() * 1e3

    def _metadata_edited(self):
        self.metadataChanged.emit(self.dt, self.umax)

    # --------------------------------------------------------------- update ---
    def set_pulse(self, pulse):
        self.pulse = pulse
        self.schedule()

    def schedule(self, *args):
        self._timer.start()

    def recompute(self):
        ps = self.pulse
        if ps is None or ps.N == 0:
            for lab in self.labels.values():
                lab.setText("–")
                lab.setStyleSheet("")
            self.labels["info"].setText("not calculated")
            self.view3d.plot_curve(None)
            self.fid_curve.setData([], [])
            return
        try:
            g = metrics.gate(ps)
            rb = metrics.robustness(ps)
            b1 = metrics.amplitude_robustness(ps)
            curve = ps_to_curve(ps)
        except Exception:
            self.status.setText("calculation failed:\n" + traceback.format_exc(limit=2))
            return
        self.labels["info"].setText(f"{ps.N} sub-pulses, {ps.duration * 1e6:.2f} µs")
        self.labels["amp"].setText(f"{np.max(ps.amplitude):.2f} % = {np.max(ps.amplitude) * ps.umax / 1e5:.3f} kHz")
        self.labels["gate"].setText(f"{g['angle_deg']:.4f}°")
        self.labels["axis"].setText(fmt_vec(g["axis"], 4))
        for k in ("A1", "A2", "A3"):
            v = rb[k + "_norm"]
            self._set_value(self.labels[k], v, f"{v:.3e}")
            self.labels[k].setToolTip(f"A{k[1]} / T^{k[1]} = {fmt_vec(rb[k] / ps.duration ** int(k[1]), 4)}")
        self._set_value(self.labels["B1"], b1["B1_norm"], f"{b1['B1_norm']:.3e}")
        self.view3d.plot_curve(curve, "error curve (start green, end red)")
        if self.fid_check.isChecked():
            half = self.fid_range.value() * 1e3 / 2
            offs = np.linspace(-half, half, 81)
            F = metrics.fidelity_vs_offset(ps, offs)
            y = np.log10(np.maximum(1 - F, 1e-16)) if self.fid_log.isChecked() else F
            self.fid_plot.setLabel("left", "log10(1-F)" if self.fid_log.isChecked() else "F")
            self.fid_curve.setData(offs / 1e3, y)
        self.last_results = {"gate": g, "robustness": rb, "B1": b1, "curve": curve}

    @staticmethod
    def _set_value(label, value, text):
        label.setText(text)
        label.setStyleSheet(f"color: {GOOD if value < ROBUST_THRESHOLD else '#e5e7eb'};")

    # --------------------------------------------------------------- output ---
    FORMATS = (".csv", ".bruker")

    def _sync_format_from_path(self, *args):
        ext = os.path.splitext(self.path_edit.text().strip())[1].lower()
        if ext in self.FORMATS:
            self.format_box.blockSignals(True)
            self.format_box.setCurrentIndex(self.FORMATS.index(ext))
            self.format_box.blockSignals(False)

    def _format_changed(self, index):
        path = self.path_edit.text().strip()
        root, ext = os.path.splitext(path)
        if not path:
            return
        self.path_edit.setText((root if ext.lower() in self.FORMATS else path) + self.FORMATS[index])

    def _browse(self):
        filters = "CSV (*.csv);;Bruker shape (*.bruker)"
        selected = "Bruker shape (*.bruker)" if self.format_box.currentIndex() == 1 else "CSV (*.csv)"
        path, chosen = QFileDialog.getSaveFileName(self, "Save pulse sequence", self.path_edit.text(), filters, selected)
        if path:
            if os.path.splitext(path)[1].lower() not in self.FORMATS:
                path += ".bruker" if "bruker" in chosen else ".csv"
            self.path_edit.setText(path)
            self._sync_format_from_path()

    def save(self):
        if self.pulse is None:
            self.status.setText("no pulse sequence to save (in Curve families: click \"Calculate now\")")
            return None
        path = self.path_edit.text().strip()
        if not path:
            self.status.setText("please enter a file path")
            return None
        if os.path.splitext(path)[1].lower() not in self.FORMATS:      # no / unknown extension: use the selected format
            path += self.FORMATS[self.format_box.currentIndex()]
            self.path_edit.setText(path)
        try:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            self.pulse.save(path)
        except Exception as e:
            self.status.setText(f"saving failed: {e}")
            return None
        self.status.setText(f"saved {self.pulse.N} sub-pulses to {path}")
        return path

    def save_and_send(self):
        path = self.save()
        if path:
            self.sendRequested.emit(path)
