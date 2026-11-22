"""Curve families tab: choose a family, edit its parameters live, get the pulse sequence."""
from __future__ import annotations

import json
import traceback

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
                               QHBoxLayout, QLabel, QLineEdit, QPushButton, QRadioButton, QScrollArea, QSlider,
                               QSpinBox, QVBoxLayout, QWidget, QSplitter)

from psa.design.conversion import curve_to_ps, ps_to_curve
from psa.design.families import FAMILIES, CustomCurve, family_by_name
from psa.editor_gui.widgets import ACCENT, BAD, MUTED_FG, Mpl3D, make_plot

SLIDER_STEPS = 1000


class ParamControl(QWidget):
    """Slider + number box for one family parameter."""
    changed = Signal()

    def __init__(self, param, parent=None):
        super().__init__(parent)
        self.param = param
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        if param.integer:
            self.slider.setRange(int(param.minimum), int(param.maximum))
            self.box = QSpinBox(); self.box.setRange(int(param.minimum), int(param.maximum))
        else:
            self.slider.setRange(0, SLIDER_STEPS)
            self.box = QDoubleSpinBox(); self.box.setRange(param.minimum, param.maximum); self.box.setDecimals(4)
            self.box.setSingleStep((param.maximum - param.minimum) / 200)
        if param.unit:
            self.box.setSuffix(" " + param.unit)
        self.box.setMinimumWidth(110)
        layout.addWidget(self.slider, stretch=1)
        layout.addWidget(self.box)
        self.slider.valueChanged.connect(self._from_slider)
        self.box.valueChanged.connect(self._from_box)
        self.set_value(param.default)

    def _to_slider(self, v):
        p = self.param
        return int(round(v)) if p.integer else int(round((v - p.minimum) / (p.maximum - p.minimum) * SLIDER_STEPS))

    def _from_slider(self, s):
        p = self.param
        v = s if p.integer else p.minimum + s / SLIDER_STEPS * (p.maximum - p.minimum)
        self.box.blockSignals(True); self.box.setValue(v); self.box.blockSignals(False)
        self.changed.emit()

    def _from_box(self, v):
        self.slider.blockSignals(True); self.slider.setValue(self._to_slider(v)); self.slider.blockSignals(False)
        self.changed.emit()

    def value(self):
        return self.box.value()

    def set_value(self, v):
        self.box.blockSignals(True); self.slider.blockSignals(True)
        self.box.setValue(v); self.slider.setValue(self._to_slider(self.box.value()))
        self.box.blockSignals(False); self.slider.blockSignals(False)


class FamiliesTab(QWidget):
    pulseChanged = Signal(object)
    editInWaveform = Signal(object)

    def __init__(self, dt=0.5e-6, umax=1e4, parent=None):
        super().__init__(parent)
        self.dt, self.umax = dt, umax
        self.pulse = None
        self.controls = {}
        outer = QHBoxLayout(self)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter)

        # ---- left: controls
        left = QWidget(); ll = QVBoxLayout(left)
        self.family_box = QComboBox()
        self.family_box.addItems([f.name for f in FAMILIES])
        ll.addWidget(self.family_box)
        self.description = QLabel(); self.description.setWordWrap(True)
        self.description.setStyleSheet(f"color: {MUTED_FG};")
        ll.addWidget(self.description)

        self.expr_group = QGroupBox("Formulas")
        ef = QFormLayout(self.expr_group)
        self.expr_edits = []
        for axis, text in zip("xyz", CustomCurve.default_expressions):
            e = QLineEdit(text)
            e.textChanged.connect(self.schedule)
            ef.addRow(f"{axis}(u) =", e)
            self.expr_edits.append(e)
        ll.addWidget(self.expr_group)

        self.param_group = QGroupBox("Parameters")
        self.param_form = QFormLayout(self.param_group)
        ll.addWidget(self.param_group)

        self.dur_group = QGroupBox("Pulse duration")
        dl = QVBoxLayout(self.dur_group)
        self.dur_short = QRadioButton("shortest possible (uses the full amplitude)")
        self.dur_fixed = QRadioButton("fixed:")
        self.dur_short.setChecked(True)
        group = QButtonGroup(self); group.addButton(self.dur_short); group.addButton(self.dur_fixed)
        self.dur_spin = QDoubleSpinBox(); self.dur_spin.setDecimals(4); self.dur_spin.setRange(1e-4, 1e4)
        self.dur_spin.setSuffix(" ms"); self.dur_spin.setValue(1.0)
        row = QHBoxLayout(); row.addWidget(self.dur_fixed); row.addWidget(self.dur_spin)
        dl.addWidget(self.dur_short); dl.addLayout(row)
        for w in (self.dur_short, self.dur_fixed):
            w.toggled.connect(self.schedule)
        self.dur_spin.valueChanged.connect(self.schedule)
        ll.addWidget(self.dur_group)

        calc = QGroupBox("Pulse calculation")
        cl = QHBoxLayout(calc)
        self.live_check = QCheckBox("calculate pulse sequence live")
        self.live_check.setChecked(True)
        self.live_check.setToolTip("Off: only the family curve is drawn while you change parameters; "
                                   "the pulse sequence and the panel values are calculated on 'Calculate now'.")
        self.calc_button = QPushButton("Calculate now")
        self.calc_button.clicked.connect(self.compute)
        self.live_check.toggled.connect(self._live_toggled)
        cl.addWidget(self.live_check); cl.addWidget(self.calc_button)
        ll.addWidget(calc)

        self.info = QLabel(""); self.info.setWordWrap(True)
        ll.addWidget(self.info)
        self.error = QLabel(""); self.error.setWordWrap(True); self.error.setStyleSheet(f"color: {BAD};")
        ll.addWidget(self.error)
        brow = QHBoxLayout()
        b_edit = QPushButton("Edit in Waveform tab"); b_edit.clicked.connect(lambda: self.pulse is not None and self.editInWaveform.emit(self.pulse))
        b_open = QPushButton("Open design…"); b_open.clicked.connect(self._open_design)
        b_save = QPushButton("Save design…"); b_save.clicked.connect(self._save_design)
        for b in (b_edit, b_open, b_save):
            brow.addWidget(b)
        ll.addLayout(brow)
        ll.addStretch()
        scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(left); scroll.setMinimumWidth(420)
        splitter.addWidget(scroll)

        # ---- right: resulting pulse
        right = QWidget(); rl = QVBoxLayout(right)
        self.curve_view = Mpl3D()
        self.curve_view.setMinimumHeight(320)
        rl.addWidget(self.curve_view, stretch=2)
        self.amp_plot = make_plot("amplitude [%]")
        self.phase_plot = make_plot("phase [°]")
        for p in (self.amp_plot, self.phase_plot):
            p.setMinimumHeight(180)
        self.phase_plot.setXLink(self.amp_plot)
        self.amp_curve = self.amp_plot.plot(stepMode="center", pen=pg.mkPen(ACCENT, width=1.5), fillLevel=0,
                                            brush=pg.mkBrush(56, 189, 248, 40))
        self.phase_curve = self.phase_plot.plot(stepMode="center", pen=pg.mkPen(ACCENT, width=1.5))
        rl.addWidget(self.amp_plot); rl.addWidget(self.phase_plot)
        splitter.addWidget(right)
        splitter.setStretchFactor(1, 2)

        self._timer = QTimer(self); self._timer.setSingleShot(True); self._timer.setInterval(150)
        self._timer.timeout.connect(self._update)
        self.family_box.currentIndexChanged.connect(self._family_changed)
        self._family_changed()

    # ---------------------------------------------------------------- state ---
    @property
    def family(self):
        return family_by_name(self.family_box.currentText())

    def values(self):
        return {name: c.value() for name, c in self.controls.items()}

    def _family_changed(self, *args):
        fam = self.family
        self.description.setText(fam.description)
        while self.param_form.rowCount():
            self.param_form.removeRow(0)
        self.controls = {}
        for p in fam.params:
            c = ParamControl(p)
            c.changed.connect(self.schedule)
            self.param_form.addRow(p.label + ":", c)
            self.controls[p.name] = c
        self.expr_group.setVisible(isinstance(fam, CustomCurve))
        self.dur_group.setVisible(fam.kind == "curve")
        self.schedule()

    def set_metadata(self, dt, umax):
        self.dt, self.umax = dt, umax
        self.schedule()

    def schedule(self, *args):
        self._timer.start()

    def _update(self):
        """Debounced reaction to any change: always redraw the family curve; calculate the pulse only live."""
        if self.live_check.isChecked():
            self.compute()
        else:
            self.draw_curve()
            self._mark_stale()

    def _live_toggled(self, on):
        if on:
            self.compute()
        else:
            self._mark_stale()

    def _mark_stale(self):
        """The shown pulse no longer matches the parameters: clear it and the panel."""
        if self.pulse is not None or not self.info.text().startswith("pulse not calculated"):
            self.pulse = None
            self.amp_curve.setData([], [])
            self.phase_curve.setData([], [])
            self.info.setText("pulse not calculated (live calculation off) – click 'Calculate now'")
            self.pulseChanged.emit(None)

    def family_points(self):
        """Points of the current curve family (raises ValueError for invalid parameters / formulas)."""
        fam = self.family
        kwargs = {"expressions": [e.text() for e in self.expr_edits]} if isinstance(fam, CustomCurve) else {}
        return fam.points(self.values(), **kwargs)

    def draw_curve(self, ps=None):
        """3D view: the family curve (curve families, own units) or the pulse's error curve (pulse families)."""
        try:
            if self.family.kind == "curve":
                self.curve_view.plot_curve(self.family_points(), "family curve (as defined, before alignment)",
                                           scale=1.0, unit="a.u.")
            else:
                ps = ps if ps is not None else self.family.pulse(self.values(), self.dt, self.umax)
                self.curve_view.plot_curve(ps_to_curve(ps), "error curve of the pulse")
        except Exception as e:
            self.error.setText(str(e))
            return False
        self.error.setText("")
        return True

    # -------------------------------------------------------------- compute ---
    def compute(self):
        fam = self.family
        try:
            if fam.kind == "curve":
                pts = self.family_points()
                duration = self.dur_spin.value() * 1e-3 if self.dur_fixed.isChecked() else None
                ps, info = curve_to_ps(pts, self.dt, self.umax, duration=duration)
                err = np.max(np.linalg.norm(ps_to_curve(ps) - info["target_points"], axis=1)) / info["duration"]
                text = (f"duration {info['duration'] * 1e6:.2f} µs ({ps.N} sub-pulses); shortest possible "
                        f"{info['min_duration'] * 1e6:.2f} µs; max. amplitude {info['max_amplitude_percent']:.2f} %; "
                        f"curve reproduced to {err:.1e} of its length")
                if not info["amplitude_ok"]:
                    self.error.setText("the duration is too short for this curve: amplitude above 100 %")
                else:
                    self.error.setText("")
            else:
                ps = fam.pulse(self.values(), self.dt, self.umax)
                text = f"duration {ps.duration * 1e6:.2f} µs ({ps.N} sub-pulses)"
                self.error.setText("")
        except Exception as e:
            self.error.setText(str(e) or traceback.format_exc(limit=1))
            return None
        self.info.setText(text)
        self.pulse = ps
        if fam.kind == "curve":
            self.curve_view.plot_curve(pts, "family curve (as defined, before alignment)", scale=1.0, unit="a.u.")
        else:
            self.curve_view.plot_curve(ps_to_curve(ps), "error curve of the pulse")
        x = np.arange(ps.N + 1) * ps.dt * 1e6
        self.amp_curve.setData(x, ps.amplitude)
        self.phase_curve.setData(x, ps.phase)
        self.pulseChanged.emit(ps)
        return ps

    def current_pulse(self):
        if self.pulse is not None or not self.live_check.isChecked():
            return self.pulse
        return self.compute()

    # -------------------------------------------------------------- designs ---
    def to_dict(self):
        return {"type": "family", "family": self.family.name, "values": self.values(),
                "expressions": [e.text() for e in self.expr_edits],
                "duration_ms": self.dur_spin.value() if self.dur_fixed.isChecked() else None}

    def load_design(self, d):
        self.family_box.setCurrentText(d["family"])
        for name, v in d.get("values", {}).items():
            if name in self.controls:
                self.controls[name].set_value(v)
        for e, text in zip(self.expr_edits, d.get("expressions", [])):
            e.setText(text)
        if d.get("duration_ms") is not None:
            self.dur_fixed.setChecked(True); self.dur_spin.setValue(d["duration_ms"])
        else:
            self.dur_short.setChecked(True)
        self._timer.stop()
        self.compute()

    def _open_design(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open family design", "", "Design (*.json)")
        if path:
            with open(path) as fh:
                self.load_design(json.load(fh))

    def _save_design(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save family design", "family_design.json", "Design (*.json)")
        if path:
            with open(path, "w") as fh:
                json.dump(self.to_dict(), fh, indent=1)
