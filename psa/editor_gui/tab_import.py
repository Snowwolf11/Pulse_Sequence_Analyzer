"""Import tab: load a curve (points or curvature / torsion) from a text file and convert it."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QButtonGroup, QComboBox, QDoubleSpinBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel,
                               QLineEdit, QPushButton, QRadioButton, QSplitter, QVBoxLayout, QWidget)

from psa.design.conversion import curve_to_ps, ps_to_curve
from psa.design.importing import FORMATS, import_curve
from psa.editor_gui.widgets import ACCENT, BAD, MUTED_FG, Mpl3D, make_plot

HELP = ("Text / CSV file, one row per line; separators: comma, semicolon, tab or spaces; header and '#' lines "
        "are skipped. Units do not matter: the curve is aligned (start at the origin, initial tangent +z) and "
        "scaled so that its length equals the pulse duration. A curvature/torsion table is integrated with the "
        "Frenet–Serret equations first.")


class ImportTab(QWidget):
    pulseChanged = Signal(object)
    editInWaveform = Signal(object)

    def __init__(self, dt=0.5e-6, umax=1e4, parent=None):
        super().__init__(parent)
        self.dt, self.umax = dt, umax
        self.points = None
        self.pulse = None
        outer = QHBoxLayout(self)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter)

        left = QWidget(); ll = QVBoxLayout(left)
        hint = QLabel(HELP); hint.setWordWrap(True); hint.setStyleSheet(f"color: {MUTED_FG};")
        ll.addWidget(hint)
        frow = QHBoxLayout()
        self.path_edit = QLineEdit(); self.path_edit.setPlaceholderText("curve file (.csv / .txt)")
        browse = QPushButton("…"); browse.setFixedWidth(30); browse.clicked.connect(self._browse)
        frow.addWidget(self.path_edit); frow.addWidget(browse)
        ll.addLayout(frow)
        self.format_box = QComboBox()
        self.format_keys = list(FORMATS)
        self.format_box.addItems([FORMATS[k] for k in self.format_keys])
        ll.addWidget(self.format_box)
        self.load_button = QPushButton("Load and convert"); self.load_button.clicked.connect(self.load)
        ll.addWidget(self.load_button)

        dur = QGroupBox("Pulse duration")
        dl = QVBoxLayout(dur)
        self.dur_short = QRadioButton("shortest possible (uses the full amplitude)"); self.dur_short.setChecked(True)
        self.dur_fixed = QRadioButton("fixed:")
        g = QButtonGroup(self); g.addButton(self.dur_short); g.addButton(self.dur_fixed)
        self.dur_spin = QDoubleSpinBox(); self.dur_spin.setDecimals(4); self.dur_spin.setRange(1e-4, 1e4)
        self.dur_spin.setSuffix(" ms"); self.dur_spin.setValue(1.0)
        row = QHBoxLayout(); row.addWidget(self.dur_fixed); row.addWidget(self.dur_spin)
        dl.addWidget(self.dur_short); dl.addLayout(row)
        self.dur_short.toggled.connect(lambda *a: self.convert())
        self.dur_spin.editingFinished.connect(self.convert)
        ll.addWidget(dur)
        self.info = QLabel(""); self.info.setWordWrap(True)
        ll.addWidget(self.info)
        self.error = QLabel(""); self.error.setWordWrap(True); self.error.setStyleSheet(f"color: {BAD};")
        ll.addWidget(self.error)
        b_edit = QPushButton("Edit in Waveform tab")
        b_edit.clicked.connect(lambda: self.pulse is not None and self.editInWaveform.emit(self.pulse))
        ll.addWidget(b_edit)
        ll.addStretch()
        splitter.addWidget(left)

        right = QWidget(); rl = QVBoxLayout(right)
        self.view3d = Mpl3D(); self.view3d.setMinimumHeight(300)
        rl.addWidget(self.view3d, stretch=2)
        self.amp_plot = make_plot("amplitude [%]"); self.phase_plot = make_plot("phase [°]")
        self.phase_plot.setXLink(self.amp_plot)
        for p in (self.amp_plot, self.phase_plot):
            p.setMinimumHeight(180)
            rl.addWidget(p)
        self.amp_curve = self.amp_plot.plot(stepMode="center", pen=pg.mkPen(ACCENT, width=1.5), fillLevel=0,
                                            brush=pg.mkBrush(56, 189, 248, 40))
        self.phase_curve = self.phase_plot.plot(stepMode="center", pen=pg.mkPen(ACCENT, width=1.5))
        splitter.addWidget(right)
        splitter.setStretchFactor(1, 2)

    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open curve", "", "Curve (*.csv *.txt *.dat);;All files (*)")
        if path:
            self.path_edit.setText(path)
            self.load()

    def load(self):
        path = self.path_edit.text().strip()
        fmt = self.format_keys[self.format_box.currentIndex()]
        try:
            self.points = import_curve(path, fmt)
        except Exception as e:
            self.error.setText(f"could not read the curve: {e}")
            self.points = None
            return None
        self.view3d.plot_curve(self.points, f"imported curve ({len(self.points)} points)", scale=1.0, unit="a.u.")
        return self.convert()

    def convert(self):
        if self.points is None:
            return None
        duration = self.dur_spin.value() * 1e-3 if self.dur_fixed.isChecked() else None
        try:
            ps, info = curve_to_ps(self.points, self.dt, self.umax, duration=duration)
        except Exception as e:
            self.error.setText(f"conversion failed: {e}")
            return None
        err = np.max(np.linalg.norm(ps_to_curve(ps) - info["target_points"], axis=1)) / info["duration"]
        self.info.setText(f"duration {info['duration'] * 1e6:.2f} µs ({ps.N} sub-pulses); shortest possible "
                          f"{info['min_duration'] * 1e6:.2f} µs; max. amplitude {info['max_amplitude_percent']:.2f} %; "
                          f"curve reproduced to {err:.1e} of its length")
        self.error.setText("" if info["amplitude_ok"] else "the duration is too short for this curve: amplitude above 100 %")
        self.pulse = ps
        x = np.arange(ps.N + 1) * ps.dt * 1e6
        self.amp_curve.setData(x, ps.amplitude)
        self.phase_curve.setData(x, ps.phase)
        self.pulseChanged.emit(ps)
        return ps

    def set_metadata(self, dt, umax):
        self.dt, self.umax = dt, umax
        self.convert()

    def current_pulse(self):
        return self.pulse
