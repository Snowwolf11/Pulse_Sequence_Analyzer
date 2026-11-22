"""Waveform tab: visual amplitude / phase editing with keyframes.

Mouse: drag a keyframe to move it (snapped to sub-pulse starts), double-click to add one,
right-click a keyframe to delete it, click to select it (then edit value / mode precisely).
"""
from __future__ import annotations

import json

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QComboBox, QDoubleSpinBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel,
                               QPushButton, QVBoxLayout, QWidget)

from psa.design.pulse import PulseSequence
from psa.design.waveform import HOLD, LINEAR, Waveform
from psa.editor_gui.widgets import ACCENT, MUTED_FG, KeyframeItem, make_plot


class WaveformTab(QWidget):
    pulseChanged = Signal(object)
    requestTextPulse = Signal()          # ask the window for the Text tab's pulse

    def __init__(self, dt=0.5e-6, umax=1e4, parent=None):
        super().__init__(parent)
        self.model = Waveform(N=int(round(100e-6 / dt)), dt=dt, umax=umax)
        self.selection = None            # (channel name, position)
        layout = QVBoxLayout(self)

        # ---- duration
        top = QHBoxLayout()
        top.addWidget(QLabel("total duration:"))
        self.duration_spin = QDoubleSpinBox(); self.duration_spin.setDecimals(4); self.duration_spin.setRange(1e-4, 1e4)
        self.duration_spin.setSuffix(" ms"); self.duration_spin.setSingleStep(0.01)
        self.duration_spin.setValue(self.model.duration * 1e3)
        self.duration_spin.editingFinished.connect(self._duration_edited)
        top.addWidget(self.duration_spin)
        self.n_label = QLabel(""); self.n_label.setStyleSheet(f"color: {MUTED_FG};")
        top.addWidget(self.n_label)
        top.addStretch()
        b_reset = QPushButton("Reset"); b_reset.clicked.connect(self.reset)
        top.addWidget(b_reset)
        layout.addLayout(top)

        # ---- import and design files
        imp = QHBoxLayout()
        b_text = QPushButton("From Text tab"); b_text.clicked.connect(self.requestTextPulse.emit)
        b_file = QPushButton("Open pulse file…"); b_file.clicked.connect(self._open_pulse_file)
        imp.addWidget(b_text); imp.addWidget(b_file)
        imp.addWidget(QLabel("tolerance:"))
        self.amp_tol = QDoubleSpinBox(); self.amp_tol.setRange(0, 100); self.amp_tol.setValue(0.5); self.amp_tol.setSuffix(" %")
        self.phase_tol = QDoubleSpinBox(); self.phase_tol.setRange(0, 360); self.phase_tol.setValue(1.0); self.phase_tol.setSuffix(" °")
        for w in (self.amp_tol, self.phase_tol):
            w.setToolTip("keyframes are placed so that the imported pulse is reproduced within this tolerance (0 = exact)")
            imp.addWidget(w)
        imp.addStretch()
        b_open = QPushButton("Open design…"); b_open.clicked.connect(self._open_design)
        b_save = QPushButton("Save design…"); b_save.clicked.connect(self._save_design)
        imp.addWidget(b_open); imp.addWidget(b_save)
        layout.addLayout(imp)

        # ---- selected keyframe
        sel = QGroupBox("Selected keyframe")
        sl = QHBoxLayout(sel)
        self.sel_label = QLabel("none (click one)")
        self.sel_value = QDoubleSpinBox(); self.sel_value.setRange(-1e7, 1e7); self.sel_value.setDecimals(4)
        self.sel_value.editingFinished.connect(self._value_edited)
        self.sel_mode = QComboBox(); self.sel_mode.addItems(["linear to next", "hold until next"])
        self.sel_mode.currentIndexChanged.connect(self._mode_edited)
        b_del = QPushButton("Delete"); b_del.clicked.connect(self._delete_selected)
        for w in (self.sel_label, QLabel("value:"), self.sel_value, QLabel("segment:"), self.sel_mode, b_del):
            sl.addWidget(w)
        sl.addStretch()
        layout.addWidget(sel)

        hint = QLabel("drag keyframes · double-click to add · right-click to delete · keyframes snap to sub-pulse starts")
        hint.setStyleSheet(f"color: {MUTED_FG};")
        layout.addWidget(hint)

        # ---- plots
        self.plots, self.steps, self.items = {}, {}, {}
        for name, label in (("amplitude", "amplitude [%]"), ("phase", "phase [°]")):
            p = make_plot(label)
            p.setMinimumHeight(220)
            step = p.plot(stepMode="center", pen=pg.mkPen(ACCENT, width=1.5), fillLevel=0 if name == "amplitude" else None,
                          brush=pg.mkBrush(56, 189, 248, 40) if name == "amplitude" else None)
            item = KeyframeItem()
            p.addItem(item)
            item.sigDragged.connect(lambda pos, x, y, n=name: self._dragged(n, pos, x, y))
            item.sigSelected.connect(lambda pos, n=name: self._select(n, pos))
            item.sigDeleteRequested.connect(lambda pos, n=name: self._delete(n, pos))
            p.scene().sigMouseClicked.connect(lambda ev, n=name: self._scene_clicked(n, ev))
            self.plots[name], self.steps[name], self.items[name] = p, step, item
            layout.addWidget(p, stretch=1)
        self.plots["phase"].setXLink(self.plots["amplitude"])
        self.plots["amplitude"].setYRange(0, 105)
        self.refresh(emit=False)

    # ------------------------------------------------------------- helpers ---
    def channel(self, name):
        return getattr(self.model, name)

    def _x_to_index(self, x_us):
        return int(round(x_us * 1e-6 / self.model.dt))

    def refresh(self, emit=True):
        m = self.model
        x_edges = np.arange(m.N + 1) * m.dt * 1e6
        for name in ("amplitude", "phase"):
            ch = self.channel(name)
            self.steps[name].setData(x_edges, ch.values(m.N))
            sel = self.selection[1] if self.selection and self.selection[0] == name else None
            xs = np.array([k.index for k in ch.keyframes]) * m.dt * 1e6
            self.items[name].set_points(xs, np.array([k.value for k in ch.keyframes]), sel)
        self.n_label.setText(f"= {m.N} sub-pulses of {m.dt * 1e6:g} µs")
        self._update_selection_widgets()
        if emit:
            self.pulseChanged.emit(self.current_pulse())

    def current_pulse(self) -> PulseSequence:
        return self.model.to_pulse()

    def set_metadata(self, dt, umax):
        self.model.umax = umax
        if abs(dt - self.model.dt) > 1e-15:
            self.model.set_dt(dt)
        self.duration_spin.setValue(self.model.duration * 1e3)
        self.selection = None
        self.refresh()

    # --------------------------------------------------------- interaction ---
    def _dragged(self, name, pos, x, y):
        self.channel(name).move(pos, self._x_to_index(x), y, self.model.N)
        self.selection = (name, pos)
        self.refresh()

    def _select(self, name, pos):
        self.selection = (name, pos)
        self.refresh(emit=False)

    def _delete(self, name, pos):
        self.channel(name).remove(pos)
        self.selection = None
        self.refresh()

    def _scene_clicked(self, name, ev):
        if not ev.double() or ev.button() != Qt.MouseButton.LeftButton:
            return
        vb = self.plots[name].getPlotItem().vb
        if not vb.sceneBoundingRect().contains(ev.scenePos()):
            return
        p = vb.mapSceneToView(ev.scenePos())
        self.add_keyframe(name, self._x_to_index(p.x()), p.y())

    def add_keyframe(self, name, index, value, mode=None):
        ch = self.channel(name)
        index = int(min(max(index, 0), self.model.N - 1))
        ch.set(index, value, mode)
        self.selection = (name, ch.find(index))
        self.refresh()

    def _update_selection_widgets(self):
        enabled = self.selection is not None
        for w in (self.sel_value, self.sel_mode):
            w.setEnabled(enabled)
        if not enabled:
            self.sel_label.setText("none (click one)")
            return
        name, pos = self.selection
        k = self.channel(name).keyframes[pos]
        self.sel_label.setText(f"{name} at {k.index * self.model.dt * 1e6:.3f} µs (sub-pulse {k.index})")
        for w, v in ((self.sel_value, k.value), (self.sel_mode, 1 if k.mode == HOLD else 0)):
            w.blockSignals(True)
            w.setValue(v) if w is self.sel_value else w.setCurrentIndex(v)
            w.blockSignals(False)

    def _value_edited(self):
        if self.selection:
            name, pos = self.selection
            k = self.channel(name).keyframes[pos]
            k.value = self.channel(name).clamp(self.sel_value.value())
            self.refresh()

    def _mode_edited(self, idx):
        if self.selection:
            name, pos = self.selection
            self.channel(name).keyframes[pos].mode = HOLD if idx == 1 else LINEAR
            self.refresh()

    def _delete_selected(self):
        if self.selection:
            self._delete(*self.selection)

    def _duration_edited(self):
        self.model.set_duration(self.duration_spin.value() * 1e-3)
        self.duration_spin.setValue(self.model.duration * 1e3)
        self.selection = None
        self.refresh()

    def reset(self):
        self.model = Waveform(self.model.N, self.model.dt, self.model.umax)
        self.selection = None
        self.refresh()

    # -------------------------------------------------------------- import ---
    def load_pulse(self, ps: PulseSequence):
        """Replace the waveform by keyframes reproducing ``ps`` (within the import tolerances).
        The editor keeps its own metadata; the pulse's sub-pulses are taken as they are."""
        self.model = Waveform.from_pulse(PulseSequence(ps.amplitude, ps.phase, self.model.dt, self.model.umax),
                                         self.amp_tol.value(), self.phase_tol.value())
        self.duration_spin.setValue(self.model.duration * 1e3)
        self.selection = None
        self.refresh()

    def _open_pulse_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open pulse sequence", "", "Pulse sequence (*.bruker *.csv *.txt);;All files (*)")
        if path:
            self.load_pulse(PulseSequence.load(path, self.model.dt, self.model.umax))

    def _open_design(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open waveform design", "", "Design (*.json)")
        if path:
            with open(path) as fh:
                self.load_design(json.load(fh))

    def load_design(self, d):
        self.model = Waveform.from_dict(d)
        self.duration_spin.setValue(self.model.duration * 1e3)
        self.selection = None
        self.refresh()

    def _save_design(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save waveform design", "waveform_design.json", "Design (*.json)")
        if path:
            with open(path, "w") as fh:
                json.dump(self.model.to_dict(), fh, indent=1)
