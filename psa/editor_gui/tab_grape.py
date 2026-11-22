"""Optimize (GRAPE) tab: optimize every sub-pulse (Cartesian controls) with exact gradients.

Start pulse: the current pulse of another tab, a rectangular pulse or a random one. Robustness: local
(Taylor orders of the detuning error, SCQC, at a list of offsets; optional first-order B1) or classical
broadband (the gate on a grid of offsets x amplitude scalings). Gate goal: none, UR or PP. The optimization
runs in a separate process (psa.design.grape.run_grape_job); the plots follow the best pulse so far.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import queue as queue_module
import time

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QGridLayout, QGroupBox,
                               QHBoxLayout, QLabel, QLineEdit, QPushButton, QRadioButton, QSpinBox, QSplitter,
                               QVBoxLayout, QWidget)

from psa.design.grape import GrapeConfig, run_grape_job, selection_key
from psa.design.optimize import SolverSettings
from psa.design.pulse import PulseSequence
from psa.editor_gui.widgets import ACCENT, ACCENT_2, BAD, GOOD, MUTED_FG, lock_heights, make_plot

START_COLORS = ("#f59e0b", "#a78bfa", "#34d399", "#f472b6", "#60a5fa", "#fbbf24", "#f87171", "#2dd4bf")
SOURCES = ("Text tab", "Waveform tab", "Curve families tab", "Curve builder tab", "Import curve tab",
           "rectangular", "random")


def vector_row(values, rng=(-1.0, 1.0)):
    row = QHBoxLayout()
    spins = []
    for name, v in zip("xyz", values):
        sp = QDoubleSpinBox(); sp.setRange(*rng); sp.setDecimals(3); sp.setSingleStep(0.1); sp.setValue(v)
        sp.setPrefix(name + " ")
        row.addWidget(sp); spins.append(sp)
    return row, spins


class GrapeTab(QWidget):
    pulseChanged = Signal(object)
    editInWaveform = Signal(object)
    requestStart = Signal(str)               # ask the window for another tab's current pulse
    settingsRequested = Signal()

    def __init__(self, dt=0.5e-6, umax=1e4, parent=None):
        super().__init__(parent)
        self.dt, self.umax = dt, umax
        self.start_pulse = None
        self.result = None
        self.settings_provider = None         # callable() -> SolverSettings
        self._process = None
        self._history = []
        outer = QHBoxLayout(self)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter)
        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(0, 0, 0, 0)

        # ---- start pulse
        sbox = QGroupBox("Start pulse")
        sl = QGridLayout(sbox)
        self.source_box = QComboBox(); self.source_box.addItems(SOURCES); self.source_box.setCurrentText("rectangular")
        self.source_box.currentTextChanged.connect(self._source_changed)
        self.load_button = QPushButton("Load start pulse"); self.load_button.clicked.connect(self.load_start)
        sl.addWidget(QLabel("from:"), 0, 0); sl.addWidget(self.source_box, 0, 1); sl.addWidget(self.load_button, 0, 2)
        self.dur_spin = QDoubleSpinBox(); self.dur_spin.setRange(0.1, 1e5); self.dur_spin.setDecimals(2)
        self.dur_spin.setSuffix(" µs"); self.dur_spin.setValue(200.0)
        self.amp_spin = QDoubleSpinBox(); self.amp_spin.setRange(0, 100); self.amp_spin.setSuffix(" %"); self.amp_spin.setValue(100.0)
        self.phase_spin = QDoubleSpinBox(); self.phase_spin.setRange(-1e5, 1e5); self.phase_spin.setSuffix(" °")
        self.seed_spin = QSpinBox(); self.seed_spin.setRange(0, 10 ** 6); self.seed_spin.setPrefix("seed ")
        self.gen_widgets = [self.dur_spin, self.amp_spin, self.phase_spin, self.seed_spin]
        sl.addWidget(QLabel("duration:"), 1, 0); sl.addWidget(self.dur_spin, 1, 1)
        sl.addWidget(self.seed_spin, 1, 2)
        sl.addWidget(QLabel("amplitude:"), 2, 0); sl.addWidget(self.amp_spin, 2, 1); sl.addWidget(self.phase_spin, 2, 2)
        self.start_info = QLabel("no start pulse loaded"); self.start_info.setStyleSheet(f"color: {MUTED_FG};")
        self.start_info.setWordWrap(True)
        sl.addWidget(self.start_info, 3, 0, 1, 3)
        ll.addWidget(sbox)

        # ---- robustness
        rbox = QGroupBox("Robustness")
        rl = QVBoxLayout(rbox)
        self.local_radio = QRadioButton("local: Taylor orders at offsets (SCQC)")
        self.bb_radio = QRadioButton("broadband: gate on an offset × B1 grid (classical GRAPE)")
        g = QButtonGroup(self); g.addButton(self.local_radio); g.addButton(self.bb_radio)
        self.local_radio.setChecked(True)
        rl.addWidget(self.local_radio)
        self.local_widget = QWidget(); lf = QFormLayout(self.local_widget); lf.setContentsMargins(18, 0, 0, 0)
        self.offsets_edit = QLineEdit("0")
        self.offsets_edit.setToolTip("comma-separated offsets in kHz; the error expansion is made around each of them")
        lf.addRow("offsets [kHz]:", self.offsets_edit)
        orow = QHBoxLayout()
        self.order_checks = {k: QCheckBox(t) for k, t in ((1, "1st"), (2, "2nd"), (3, "3rd"))}
        self.order_checks[1].setChecked(True)
        for c in self.order_checks.values():
            orow.addWidget(c)
        self.b1_check = QCheckBox("B1 (1st)")
        self.b1_check.setToolTip("first-order robustness against amplitude scaling (at every offset)")
        orow.addWidget(self.b1_check); orow.addStretch()
        lf.addRow("orders:", orow)
        rl.addWidget(self.local_widget)
        rl.addWidget(self.bb_radio)
        self.bb_widget = QWidget(); bg = QFormLayout(self.bb_widget); bg.setContentsMargins(18, 0, 0, 0)
        self.bb_spins = {}
        for key, label, lo, hi, n, unit, rng in (("offset", "offset", -5.0, 5.0, 11, " kHz", (-1e4, 1e4)),
                                                 ("scale", "B1 scale", 1.0, 1.0, 1, "", (0.0, 10.0))):
            a, b = QDoubleSpinBox(), QDoubleSpinBox()
            for sp, v in ((a, lo), (b, hi)):
                sp.setRange(*rng); sp.setDecimals(3); sp.setSuffix(unit); sp.setValue(v)
            c = QSpinBox(); c.setRange(1, 1000); c.setValue(n); c.setSuffix(" points")
            row = QHBoxLayout(); row.addWidget(a); row.addWidget(QLabel("to")); row.addWidget(b)
            bg.addRow(label + ":", row)
            bg.addRow("", c)
            self.bb_spins[key] = (a, b, c)
        rl.addWidget(self.bb_widget)
        self.local_radio.toggled.connect(self._mode_changed)
        ll.addWidget(rbox)

        # ---- gate
        gbox = QGroupBox("Gate goal")
        gl = QVBoxLayout(gbox)
        self.gate_box = QComboBox(); self.gate_box.addItems(["none", "UR (universal rotation)", "PP (point to point)"])
        self.gate_box.setCurrentIndex(1)
        self.gate_box.currentIndexChanged.connect(self._gate_changed)
        gl.addWidget(self.gate_box)
        self.ur_widget = QWidget(); uf = QFormLayout(self.ur_widget); uf.setContentsMargins(0, 0, 0, 0)
        self.angle_spin = QDoubleSpinBox(); self.angle_spin.setRange(-720, 720); self.angle_spin.setDecimals(3)
        self.angle_spin.setSuffix(" °"); self.angle_spin.setValue(180.0)
        uf.addRow("angle:", self.angle_spin)
        row, self.axis_spins = vector_row((1.0, 0.0, 0.0)); uf.addRow("axis:", row)
        gl.addWidget(self.ur_widget)
        self.pp_widget = QWidget(); pf = QFormLayout(self.pp_widget); pf.setContentsMargins(0, 0, 0, 0)
        row, self.m0_spins = vector_row((0.0, 0.0, 1.0)); pf.addRow("from m0:", row)
        row, self.mt_spins = vector_row((1.0, 0.0, 0.0)); pf.addRow("to:", row)
        gl.addWidget(self.pp_widget)
        ll.addWidget(gbox)

        # ---- controls and smoothness
        cbox = QGroupBox("Controls and smoothness")
        cf = QFormLayout(cbox)
        crow = QHBoxLayout()
        self.full_radio = QRadioButton("every sub-pulse"); self.spline_radio = QRadioButton("smooth spline,")
        cg = QButtonGroup(self); cg.addButton(self.full_radio); cg.addButton(self.spline_radio)
        self.full_radio.setChecked(True)
        self.full_radio.setToolTip("classical GRAPE: u_x, u_y of every sub-pulse are free")
        self.spline_radio.setToolTip("u_x, u_y are cubic B-splines: smooth by construction, fewer variables")
        self.spline_spin = QSpinBox(); self.spline_spin.setRange(4, 100000); self.spline_spin.setValue(20)
        self.spline_spin.setSuffix(" control points")
        crow.addWidget(self.spline_radio); crow.addWidget(self.spline_spin); crow.addStretch()
        cf.addRow("optimize:", self.full_radio)
        cf.addRow("", crow)
        self.spline_radio.toggled.connect(lambda on: self.spline_spin.setEnabled(on))
        self.spline_spin.setEnabled(False)
        self.smooth_edit = QLineEdit("0")
        self.smooth_edit.setToolTip("weight w of the smoothness term w·N·Σ|u(k+1) − u(k)|² (soft; lowered step by step "
                                    "until the goals are reached exactly). Try 1e-4. 0 = off")
        cf.addRow("smoothness weight:", self.smooth_edit)
        self.amp_rate_spin = QDoubleSpinBox(); self.amp_rate_spin.setRange(0, 1e6); self.amp_rate_spin.setDecimals(2)
        self.amp_rate_spin.setSuffix(" %/µs"); self.amp_rate_spin.setSpecialValueText("off")
        self.amp_rate_spin.setToolTip("soft limit of the amplitude change rate (no cost below it); off = 0")
        self.phase_rate_spin = QDoubleSpinBox(); self.phase_rate_spin.setRange(0, 1e6); self.phase_rate_spin.setDecimals(2)
        self.phase_rate_spin.setSuffix(" kHz"); self.phase_rate_spin.setSpecialValueText("off")
        self.phase_rate_spin.setToolTip("soft limit of the phase change rate dφ/dt / 2π (= torsion of the error curve / 2π), "
                                        "weighted by the amplitude: phase jumps at zero amplitude are free; off = 0")
        cf.addRow("max. amplitude rate:", self.amp_rate_spin)
        cf.addRow("max. phase rate:", self.phase_rate_spin)
        reg = QHBoxLayout()
        self.reg_box = QComboBox(); self.reg_box.addItems(["none", "stay close to start"])
        self.reg_weight = QLineEdit("1e-3"); self.reg_weight.setToolTip("weight λ of λ·|u − u_start|²")
        reg.addWidget(QLabel("weight:")); reg.addWidget(self.reg_weight)
        cf.addRow("regularization:", self.reg_box)
        cf.addRow("", reg)
        ll.addWidget(cbox)

        # ---- exploration and solver
        obox = QGroupBox("Exploration and solver")
        of = QFormLayout(obox)
        mrow = QHBoxLayout()
        self.starts_spin = QSpinBox(); self.starts_spin.setRange(1, 64); self.starts_spin.setValue(1)
        self.starts_spin.setToolTip("several optimizations from perturbed starts, in parallel processes; the best "
                                    "result is kept (reached goals first, then the smoothest / lowest residual)")
        self.spread_spin = QDoubleSpinBox(); self.spread_spin.setRange(0, 10); self.spread_spin.setDecimals(3)
        self.spread_spin.setSingleStep(0.1); self.spread_spin.setValue(0.5)
        self.spread_spin.setToolTip("random perturbation of starts 2, 3, … (in the optimizer's variables; 0.5 is a "
                                    "large change, 0.05 a small one)")
        mrow.addWidget(self.starts_spin); mrow.addWidget(QLabel("spread:")); mrow.addWidget(self.spread_spin); mrow.addStretch()
        of.addRow("starts:", mrow)
        self.noise_spin = QDoubleSpinBox(); self.noise_spin.setRange(0, 1); self.noise_spin.setDecimals(4)
        self.noise_spin.setSingleStep(0.005); self.noise_spin.setValue(0.01)
        self.noise_spin.setToolTip("small seeded random change of the first start (breaks the symmetry of e.g. "
                                   "constant-phase pulses, which otherwise can trap the optimizer); 0 = off")
        of.addRow("start noise:", self.noise_spin)
        brow = QHBoxLayout()
        self.solve_button = QPushButton("Optimize"); self.solve_button.clicked.connect(self.toggle)
        self.settings_button = QPushButton("Settings…")
        self.settings_button.setToolTip("shared solver settings: evaluation limit (per start), time limit, tolerances, "
                                        "goal weights")
        self.settings_button.clicked.connect(self.settingsRequested.emit)
        brow.addWidget(self.solve_button); brow.addWidget(self.settings_button); brow.addStretch()
        of.addRow(brow)
        ll.addWidget(obox)
        self.status = QLabel(""); self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        ll.addWidget(self.status)
        arow = QHBoxLayout()
        self.restart_button = QPushButton("Result → start"); self.restart_button.clicked.connect(self.result_as_start)
        self.restart_button.setToolTip("continue from the result (e.g. with other goals)")
        self.edit_button = QPushButton("Edit in Waveform tab")
        self.edit_button.clicked.connect(lambda: self.current_pulse() is not None and self.editInWaveform.emit(self.current_pulse()))
        arow.addWidget(self.restart_button); arow.addWidget(self.edit_button); arow.addStretch()
        ll.addLayout(arow)
        ll.addStretch()
        splitter.addWidget(left)

        # ---- plots
        right = QWidget(); pl = QVBoxLayout(right); pl.setContentsMargins(0, 0, 0, 0)
        self.amp_plot, self.phase_plot = make_plot("amplitude [%]"), make_plot("phase [°]")
        self.phase_plot.setXLink(self.amp_plot)
        self.conv_plot = make_plot("residual", x_label="evaluations")
        self.conv_plot.setLogMode(y=True)
        start_pen = pg.mkPen(MUTED_FG, width=1, style=Qt.PenStyle.DashLine)
        self.amp_start = self.amp_plot.plot(stepMode="center", pen=start_pen)
        self.phase_start = self.phase_plot.plot(stepMode="center", pen=start_pen)
        self.amp_curve = self.amp_plot.plot(stepMode="center", pen=pg.mkPen(ACCENT, width=1.5), fillLevel=0,
                                            brush=pg.mkBrush(56, 189, 248, 40))
        self.phase_curve = self.phase_plot.plot(stepMode="center", pen=pg.mkPen(ACCENT, width=1.5))
        self.conv_curves = {}
        self.amp_plot.setYRange(0, 105)
        for p, h in ((self.amp_plot, 200), (self.phase_plot, 200), (self.conv_plot, 150)):
            p.setMinimumHeight(h); p.setMinimumWidth(280)
            pl.addWidget(p, stretch=2 if h == 200 else 1)
        legend = QLabel("dashed: start pulse · blue: best pulse so far / result · residual: one line per start")
        legend.setStyleSheet(f"color: {MUTED_FG};")
        legend.setWordWrap(True)
        pl.addWidget(legend)
        splitter.addWidget(right)
        splitter.setStretchFactor(1, 2)

        self._mode_changed(); self._gate_changed(); self._source_changed(self.source_box.currentText())
        lock_heights(self)
        self.load_start()

    # ------------------------------------------------------------- inputs ---
    def _mode_changed(self, *args):
        local = self.local_radio.isChecked()
        self.local_widget.setEnabled(local)
        self.bb_widget.setEnabled(not local)

    def _gate_changed(self, *args):
        i = self.gate_box.currentIndex()
        self.ur_widget.setVisible(i == 1)
        self.pp_widget.setVisible(i == 2)

    def _source_changed(self, text):
        gen = text in ("rectangular", "random")
        for w in self.gen_widgets:
            w.setEnabled(gen)
        self.amp_spin.setEnabled(text == "rectangular")
        self.phase_spin.setEnabled(text == "rectangular")
        self.seed_spin.setEnabled(text == "random")

    def config(self) -> GrapeConfig:
        """GrapeConfig from the inputs (raises ValueError with a readable message)."""
        try:
            offsets = [float(t) * 1e3 for t in self.offsets_edit.text().replace(";", ",").split(",") if t.strip()]
        except ValueError:
            raise ValueError(f"offsets: '{self.offsets_edit.text()}' is not a comma-separated list of numbers")
        try:
            weight = float(self.reg_weight.text().replace(",", "."))
        except ValueError:
            raise ValueError(f"regularization weight: '{self.reg_weight.text()}' is not a number")
        if weight < 0:
            raise ValueError("the regularization weight must not be negative")
        try:
            smooth = float(self.smooth_edit.text().replace(",", "."))
        except ValueError:
            raise ValueError(f"smoothness weight: '{self.smooth_edit.text()}' is not a number")
        if smooth < 0:
            raise ValueError("the smoothness weight must not be negative")
        (olo, ohi, on), (slo, shi, sn) = self.bb_spins["offset"], self.bb_spins["scale"]
        cfg = GrapeConfig(
            mode="local" if self.local_radio.isChecked() else "broadband", offsets_hz=offsets,
            orders=[k for k, c in self.order_checks.items() if c.isChecked()], b1=self.b1_check.isChecked(),
            bb_offset_min_hz=olo.value() * 1e3, bb_offset_max_hz=ohi.value() * 1e3, bb_offset_n=on.value(),
            bb_scale_min=slo.value(), bb_scale_max=shi.value(), bb_scale_n=sn.value(),
            gate=("none", "UR", "PP")[self.gate_box.currentIndex()],
            target_axis=[s.value() for s in self.axis_spins], target_angle_deg=self.angle_spin.value(),
            m0=[s.value() for s in self.m0_spins], m_target=[s.value() for s in self.mt_spins],
            regularization="stay_close" if self.reg_box.currentIndex() == 1 else "none", regularization_weight=weight,
            start_noise=self.noise_spin.value(), controls="spline" if self.spline_radio.isChecked() else "full",
            n_spline=self.spline_spin.value(), smooth_weight=smooth, max_amp_rate=self.amp_rate_spin.value(),
            max_phase_rate_khz=self.phase_rate_spin.value(), n_starts=self.starts_spin.value(),
            multistart_noise=self.spread_spin.value())
        if cfg.controls == "spline" and self.start_pulse is not None and cfg.n_spline > self.start_pulse.N:
            raise ValueError(f"more spline control points ({cfg.n_spline}) than sub-pulses ({self.start_pulse.N})")
        cfg.check()
        return cfg

    # -------------------------------------------------------------- start ---
    def load_start(self):
        src = self.source_box.currentText()
        if src in ("rectangular", "random"):
            n = max(1, int(round(self.dur_spin.value() * 1e-6 / self.dt)))
            if src == "rectangular":
                ps = PulseSequence(np.full(n, self.amp_spin.value()), np.full(n, self.phase_spin.value()), self.dt, self.umax)
            else:
                rng = np.random.default_rng(self.seed_spin.value())
                ps = PulseSequence(rng.uniform(0, 100, n), rng.uniform(-180, 180, n), self.dt, self.umax)
            self.set_start(ps)
        else:
            self.requestStart.emit(src)          # the window answers with set_start (or start_failed)

    def start_failed(self, text):
        self._show(text, BAD)

    def set_start(self, ps):
        if ps is None or ps.N == 0:
            self.start_failed("that tab has no valid pulse")
            return
        self.start_pulse = PulseSequence(np.asarray(ps.amplitude, float), np.asarray(ps.phase, float), ps.dt, ps.umax)
        self.result = None
        self._clear_convergence()
        self.start_info.setText(f"{ps.N} sub-pulses, {ps.duration * 1e6:.2f} µs, max. amplitude {np.max(ps.amplitude):.1f} %"
                                + ("  (amplitudes above 100 % are clipped)" if np.max(ps.amplitude) > 100 else ""))
        x = self._times(ps)
        self.amp_start.setData(x, ps.amplitude); self.phase_start.setData(x, ps.phase)
        self.amp_curve.setData([], []); self.phase_curve.setData([], [])
        self._show("")
        self.pulseChanged.emit(self.start_pulse)

    def result_as_start(self):
        if self.result is not None:
            self.set_start(self.result)

    @staticmethod
    def _times(ps):
        return np.arange(ps.N + 1) * ps.dt * 1e6

    def _clear_convergence(self):
        self._history = {}
        for c in self.conv_curves.values():
            self.conv_plot.removeItem(c)
        self.conv_curves = {}

    def _conv_curve(self, i):
        if i not in self.conv_curves:
            self.conv_curves[i] = self.conv_plot.plot(pen=pg.mkPen(START_COLORS[i % len(START_COLORS)], width=1.5))
        return self.conv_curves[i]

    def _show_result(self, ps):
        x = self._times(ps)
        self.amp_curve.setData(x, ps.amplitude); self.phase_curve.setData(x, ps.phase)

    def current_pulse(self):
        return self.result if self.result is not None else self.start_pulse

    def set_metadata(self, dt, umax):
        self.dt, self.umax = dt, umax
        if self.start_pulse is not None and not self.solving:
            ps = self.current_pulse()
            self.set_start(PulseSequence(ps.amplitude, ps.phase, dt, umax))

    # ------------------------------------------------------------ process ---
    @property
    def solving(self):
        return self._process is not None

    def toggle(self):
        if self.solving:
            self.stop()
        else:
            self.start()

    def _inputs(self):
        return [self.source_box, self.load_button, self.restart_button, self.settings_button, self.local_radio,
                self.bb_radio, self.local_widget, self.bb_widget, self.gate_box, self.ur_widget, self.pp_widget,
                self.reg_box, self.reg_weight, self.noise_spin, self.full_radio, self.spline_radio, self.spline_spin,
                self.smooth_edit, self.amp_rate_spin, self.phase_rate_spin, self.starts_spin, self.spread_spin] + self.gen_widgets

    def _set_inputs_enabled(self, on):
        for w in self._inputs():
            w.setEnabled(on)
        if on:
            self._mode_changed(); self._source_changed(self.source_box.currentText())
            self.spline_spin.setEnabled(self.spline_radio.isChecked())

    @staticmethod
    def max_parallel():
        return max(1, (os.cpu_count() or 2) - 1)

    def start(self):
        if self.start_pulse is None:
            self._show("load a start pulse first", BAD)
            return
        try:
            cfg = self.config()
        except ValueError as e:
            self._show(str(e), BAD)
            return
        settings = (self.settings_provider() if self.settings_provider else None) or SolverSettings()
        ps = self.start_pulse
        self._cfg = cfg
        self._job = {"pulse": [np.minimum(ps.amplitude, 100.0).tolist(), ps.phase.tolist(), ps.dt, ps.umax],
                     "config": cfg.to_dict(), "settings": settings.to_dict()}
        self._limit, self._started = settings.max_nfev, time.time()
        self._clear_convergence()
        self._pending = list(range(int(cfg.n_starts)))
        self._running, self._done, self._errors, self._current = {}, {}, {}, {}
        ctx = mp.get_context("spawn")
        self._ctx = ctx
        self._queue, self._stop_event = ctx.Queue(), ctx.Event()
        self._process = True
        self._launch()
        self._stop_requested_at = None
        self._poll = QTimer(self); self._poll.setInterval(100); self._poll.timeout.connect(self._poll_queue)
        self._poll.start()
        self.solve_button.setText("Stop")
        self._set_inputs_enabled(False)
        self._show("starting optimizer …")

    def _launch(self):
        while self._pending and len(self._running) < self.max_parallel() and not self._stop_event.is_set():
            i = self._pending.pop(0)
            job = dict(self._job, start_index=i)
            proc = self._ctx.Process(target=run_grape_job, args=(job, self._queue, self._stop_event), daemon=True)
            proc.start()
            self._running[i] = proc

    def stop(self, wait=False):
        if not self.solving:
            return
        self._stop_event.set()
        self._pending = []
        if self._stop_requested_at is None:
            self._stop_requested_at = time.time()
        self._show("stopping …")
        if wait:
            for proc in self._running.values():
                proc.join(5)
                if proc.is_alive():
                    proc.terminate()
            self._finish()

    def _progress_text(self):
        n = int(self._cfg.n_starts)
        best = min(self._current.values()) if self._current else float("nan")
        evals = sum(self._history.get(i, ([0], None))[0][-1] for i in self._history) if self._history else 0
        multi = f"start {len(self._done) + len(self._errors)} done / {n}, {len(self._running)} running, " if n > 1 else ""
        return (f"optimizing … {multi}{evals} evaluations (limit {self._limit} per start), "
                f"{time.time() - self._started:.0f} s, best goal residual {best:.2e}  (click Stop to interrupt)")

    def _poll_queue(self):
        changed = False
        while True:
            try:
                msg = self._queue.get_nowait()
            except queue_module.Empty:
                break
            except (EOFError, OSError):
                msg = ("error", -1, "an optimizer process ended unexpectedly")
            kind, i = msg[0], msg[1]
            if kind == "progress":
                n, c = self._history.setdefault(i, ([], []))
                n.append(msg[2]); c.append(msg[3])
                self._conv_curve(i).setData(np.array(n, float), np.maximum(np.array(c), 1e-300))
                self._current[i] = msg[3]
                changed = True
            elif kind == "snapshot":
                self._current[i] = msg[4]
                if msg[4] <= min(self._current.values()):           # show the start that is best right now
                    ps = PulseSequence(np.asarray(msg[2]), np.asarray(msg[3]), self.start_pulse.dt, self.start_pulse.umax)
                    self.result = ps
                    self._show_result(ps)
                    self.pulseChanged.emit(ps)
            elif kind == "done":
                info = msg[4]
                self._done[i] = (msg[2], msg[3], info)
                self._current[i] = info["goal_residual"]
                hist = np.asarray(info.get("history", []), float)
                if len(hist):
                    offset = max(0, info["nfev"] - len(hist))
                    self._conv_curve(i).setData(np.arange(1, len(hist) + 1) + offset, np.maximum(hist, 1e-300))
                    self._history[i] = ([info["nfev"]], [info["goal_residual"]])
                self._reap(i)
                changed = True
            elif kind == "error":
                self._errors[i] = msg[2]
                self._reap(i)
                changed = True
        if self._stop_requested_at is not None and time.time() - self._stop_requested_at > 5:
            for proc in self._running.values():
                proc.terminate()
            self._running = {}
        for i, proc in list(self._running.items()):          # crashed without a message
            if not proc.is_alive() and self._queue.empty():
                proc.join(0.1)
                if i in self._running:
                    self._errors[i] = "the optimizer process ended without a result"
                    self._reap(i)
        if not self._running and not self._pending:
            self._complete()
        elif changed:
            self._show(self._progress_text())

    def _reap(self, i):
        proc = self._running.pop(i, None)
        if proc is not None:
            proc.join(1)
        self._launch()

    def _complete(self):
        self._finish()
        if not self._done:
            err = next(iter(self._errors.values()), "optimizer terminated (no result)")
            self._show(err, BAD)
            return
        i, (amp, phase, info) = min(self._done.items(), key=lambda kv: selection_key(self._cfg, kv[1][2]))
        self.result = PulseSequence(np.asarray(amp), np.asarray(phase), self.start_pulse.dt, self.start_pulse.umax)
        self._show_result(self.result)
        for j, curve in self.conv_curves.items():
            curve.setPen(pg.mkPen(START_COLORS[j % len(START_COLORS)], width=3 if j == i else 1))
        self.last_infos = {j: v[2] for j, v in self._done.items()}
        self._report(info, best_index=i)
        self.pulseChanged.emit(self.result)

    def _finish(self):
        if getattr(self, "_poll", None) is not None:
            self._poll.stop()
            self._poll = None
        self._process = None
        self.solve_button.setText("Optimize")
        self._set_inputs_enabled(True)

    def _report(self, info, best_index=0):
        reason = info.get("stop_reason")
        if info["reached"]:
            state = "all goals reached"
        else:
            state = {"user": "stopped – best so far", "limit": "evaluation limit reached – best so far",
                     "time": "time limit reached – best so far"}.get(reason, "converged, goals not fully reached")
        fmt = lambda d: ", ".join(f"{k} = {v:.3g}" for k, v in d.items())
        n = int(self._cfg.n_starts)
        if n > 1:
            solved = sum(1 for v in self._done.values() if v[2]["reached"])
            state += f" – best of {n} starts: start {best_index + 1} ({solved} reached the goals)"
        self._show(f"{state} ({info['nfev']} evaluations, {info.get('elapsed_s', 0):.1f} s)\nresult: {fmt(info['values'])}"
                   f"\nstart: {fmt(info.get('start_summary', {}))}", GOOD if info["reached"] else BAD)

    def _show(self, text, color=None):
        self.status.setStyleSheet(f"color: {color};" if color else "")
        self.status.setText(text)
        self._fit_status()

    def _fit_status(self):
        # a word-wrapped label under-reports its height: reserve the wrapped height so that the layout
        # (inside the scroll area) grows instead of squeezing the other inputs
        width = self.status.width() if self.status.width() > 50 else 400
        self.status.setMinimumHeight(self.status.heightForWidth(width) if self.status.text() else 0)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._fit_status()
