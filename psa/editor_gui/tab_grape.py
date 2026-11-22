"""Optimize (GRAPE) tab: optimize every sub-pulse (Cartesian controls) with exact gradients.

Start pulse: the current pulse of another tab, a rectangular pulse or a random one. Robustness: local
(Taylor orders of the detuning error, SCQC, at a list of offsets; optional first-order B1) or classical
broadband (the gate on a grid of offsets x amplitude scalings). Gate goal: none, UR or PP. Exploration: a
preset (exploit / balanced / explore) or custom values (Advanced…) of the population search
(psa.design.grape.Search); its tasks run in worker processes (grape_worker), the plots follow the best pulse.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import queue as queue_module
import time

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
                               QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QPushButton, QRadioButton,
                               QSpinBox, QSplitter, QVBoxLayout, QWidget)

from psa.design.grape import EXPLORATION_PRESETS, GrapeConfig, Search, exploration_preset, grape_worker
from psa.design.optimize import SolverSettings
from psa.design.pulse import PulseSequence
from psa.editor_gui.widgets import ACCENT, ACCENT_2, BAD, GOOD, MUTED_FG, lock_heights, make_plot

START_COLORS = ("#f59e0b", "#a78bfa", "#34d399", "#f472b6", "#60a5fa", "#fbbf24", "#f87171", "#2dd4bf")
SOURCES = ("Text tab", "Waveform tab", "Curve families tab", "Curve builder tab", "Import curve tab",
           "rectangular", "random")


PRESETS = ("exploit", "balanced", "explore", "custom")
PRESET_TIPS = {
    "exploit": "one start, refined to the end (fast; finds the nearest solution)",
    "balanced": "8 smooth perturbations of the start, 2 screening rounds (8 → 4 → 2), then 4 kicks of the best",
    "explore": "24 starts (half random, smooth spread), coarse-to-fine (8, 16 spline points), 3 screening rounds "
               "(24 → 12 → 6 → 3), then 8 kicks with temperature",
    "custom": "your own values (Advanced…)",
}
DEAD_PEN = pg.mkPen("#6b7280", width=0.7)


class ExplorationDialog(QDialog):
    """Advanced exploration settings of the GRAPE tab (non-modal; values are read when Optimize starts)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("GRAPE exploration")
        f = QFormLayout(self)
        intro = QLabel("Starts: how different they are. Screening: many cheap starts, the best fraction continues with "
                       "twice the budget each round (more starts / smaller first budget = more exploration). Kicks: "
                       "perturb the best result and re-optimize. Coarse to fine: first only the large-scale shape.")
        intro.setWordWrap(True); intro.setStyleSheet(f"color: {MUTED_FG};")
        f.addRow(intro)
        self.pool_box = QComboBox(); self.pool_box.addItems(["perturb the start", "random pulses", "mix"])
        self.pool_box.setToolTip("starts 2, 3, …: perturbations of the start pulse (by the spread), random pulses, or "
                                 "alternately both")
        f.addRow("starts 2, 3, … :", self.pool_box)
        self.corr_spin = QDoubleSpinBox(); self.corr_spin.setRange(0, 1e5); self.corr_spin.setDecimals(2)
        self.corr_spin.setSuffix(" µs"); self.corr_spin.setSpecialValueText("white")
        self.corr_spin.setToolTip("correlation time of the start and kick noise: smooth noise changes the large-scale "
                                  "shape, white noise mostly adds spikes; white = 0")
        f.addRow("noise correlation:", self.corr_spin)
        self.rounds_spin = QSpinBox(); self.rounds_spin.setRange(0, 20); self.rounds_spin.setSpecialValueText("off")
        self.rounds_spin.setToolTip("successive halving: after each round only the best fraction continues; off = every "
                                    "start runs to the end")
        f.addRow("screening rounds:", self.rounds_spin)
        self.screen_spin = QSpinBox(); self.screen_spin.setRange(2, 10 ** 6); self.screen_spin.setSuffix(" evaluations")
        self.screen_spin.setToolTip("budget per candidate in round 1 (then ÷ kept fraction per round)")
        f.addRow("first round:", self.screen_spin)
        self.keep_spin = QDoubleSpinBox(); self.keep_spin.setRange(0.05, 0.95); self.keep_spin.setSingleStep(0.05)
        self.keep_spin.setDecimals(2)
        self.keep_spin.setToolTip("fraction of the candidates kept after each round")
        f.addRow("kept fraction:", self.keep_spin)
        self.cross_spin = QDoubleSpinBox(); self.cross_spin.setRange(0, 1); self.cross_spin.setSingleStep(0.1)
        self.cross_spin.setDecimals(2); self.cross_spin.setSpecialValueText("off")
        self.cross_spin.setToolTip("experimental (genetic): this fraction of the dropped places is refilled with children "
                                   "of two survivors (a blend, or the first part of one and the rest of the other); off = 0")
        f.addRow("crossover (experimental):", self.cross_spin)
        self.kicks_spin = QSpinBox(); self.kicks_spin.setRange(0, 1000); self.kicks_spin.setSpecialValueText("off")
        self.kicks_spin.setToolTip("after the starts: perturb the best result and re-optimize, keep it if better "
                                   "(skipped once the goals are reached, unless smoothing is on: then smoother wins)")
        f.addRow("kicks:", self.kicks_spin)
        self.kick_size_spin = QDoubleSpinBox(); self.kick_size_spin.setRange(0, 10); self.kick_size_spin.setDecimals(3)
        self.kick_size_spin.setSingleStep(0.05)
        self.kick_size_spin.setToolTip("size of a kick (same units as the spread)")
        f.addRow("kick size:", self.kick_size_spin)
        self.temp_spin = QDoubleSpinBox(); self.temp_spin.setRange(0, 10); self.temp_spin.setDecimals(2)
        self.temp_spin.setSingleStep(0.1); self.temp_spin.setSpecialValueText("off")
        self.temp_spin.setToolTip("accept a worse kick with probability exp(−Δ / T), Δ = increase of log10(residual) "
                                  "(or of the roughness); lets the kicks leave a local minimum; off = only improvements")
        f.addRow("kick temperature:", self.temp_spin)
        self.coarse_edit = QLineEdit("")
        self.coarse_edit.setToolTip("comma-separated spline control points of coarse stages run first (e.g. 8, 16); "
                                    "empty = off")
        f.addRow("coarse to fine:", self.coarse_edit)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.hide)
        f.addRow(buttons)
        self.widgets = [self.pool_box, self.corr_spin, self.rounds_spin, self.screen_spin, self.keep_spin,
                        self.cross_spin, self.kicks_spin, self.kick_size_spin, self.temp_spin, self.coarse_edit]

    def coarse_levels(self):
        try:
            return [int(t) for t in self.coarse_edit.text().replace(";", ",").split(",") if t.strip()]
        except ValueError:
            raise ValueError(f"coarse to fine: '{self.coarse_edit.text()}' is not a comma-separated list of integers")


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
        self._workers = []
        self._history = {}
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
        lrow = QHBoxLayout()
        self.reweight_spin = QSpinBox(); self.reweight_spin.setRange(0, 1000); self.reweight_spin.setSpecialValueText("off")
        self.reweight_spin.setSuffix(" rounds")
        self.reweight_spin.setToolTip("Lawson reweighting (both modes): if the goals are not reached at every point, each "
                                      "round raises the weights of the points (offsets / grid points) with the largest "
                                      "errors and re-optimizes. Evens out the errors: a smaller worst case at the cost "
                                      "of the mean. The round with the smallest worst-point error is kept; off = 0")
        self.reweight_power_spin = QDoubleSpinBox(); self.reweight_power_spin.setRange(0.05, 4); self.reweight_power_spin.setDecimals(2)
        self.reweight_power_spin.setSingleStep(0.1); self.reweight_power_spin.setValue(0.5)
        self.reweight_power_spin.setToolTip("strength p of the update w ← w·(error / mean error)^p: 1 is classic Lawson "
                                            "(fast, may oscillate), 0.5 is damped")
        lrow.addWidget(QLabel("even out points:")); lrow.addWidget(self.reweight_spin)
        lrow.addWidget(QLabel("strength:")); lrow.addWidget(self.reweight_power_spin); lrow.addStretch()
        rl.addLayout(lrow)
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
                                    "until the goals are reached exactly). Try 1e-4; if no start reaches the goals, 1e-5 (short pulses need faster "
                                    "modulation). 0 = off")
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
        self.limits_late_check = QCheckBox("limits only in a last stage (recommended)")
        self.limits_late_check.setChecked(True)
        self.limits_late_check.setToolTip("optimize freely first, then switch the rate limits on: limits from the start "
                                          "can block the way to good solutions (which often pass through rough "
                                          "intermediate pulses); the result keeps the limits either way")
        cf.addRow("", self.limits_late_check)
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
        self.explore_dialog = ExplorationDialog(self)
        prow = QHBoxLayout()
        self.preset_box = QComboBox(); self.preset_box.addItems(PRESETS)
        for i, name in enumerate(PRESETS):
            self.preset_box.setItemData(i, PRESET_TIPS[name], Qt.ItemDataRole.ToolTipRole)
        self.preset_box.setToolTip("exploitation ↔ exploration; the details are under Advanced…")
        self.preset_box.currentTextChanged.connect(self._preset_changed)
        self.advanced_button = QPushButton("Advanced…"); self.advanced_button.clicked.connect(self._open_advanced)
        prow.addWidget(self.preset_box); prow.addWidget(self.advanced_button); prow.addStretch()
        of.addRow("exploration:", prow)
        mrow = QHBoxLayout()
        self.starts_spin = QSpinBox(); self.starts_spin.setRange(1, 64); self.starts_spin.setValue(1)
        self.starts_spin.setRange(1, 1000)
        self.starts_spin.setToolTip("several optimizations from different starts, in parallel processes; the best "
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

        self._applying = False
        for wdg in [self.starts_spin, self.spread_spin] + self.explore_dialog.widgets:
            sig = getattr(wdg, "valueChanged", None) or getattr(wdg, "currentIndexChanged", None) or wdg.textEdited
            sig.connect(self._mark_custom)
        self._mode_changed(); self._gate_changed(); self._source_changed(self.source_box.currentText())
        self.apply_preset("exploit")
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

    # -------------------------------------------------------- exploration ---
    def _duration_s(self):
        return self.start_pulse.duration if self.start_pulse is not None else self.dur_spin.value() * 1e-6

    def apply_preset(self, name):
        """Set the exploration inputs to a preset (the noise correlation scales with the pulse duration)."""
        if name not in EXPLORATION_PRESETS:
            return
        d = exploration_preset(name, self._duration_s())
        dlg = self.explore_dialog
        self._applying = True
        try:
            self.preset_box.setCurrentText(name)
            self.starts_spin.setValue(d["n_starts"]); self.spread_spin.setValue(d["multistart_noise"])
            dlg.pool_box.setCurrentIndex(("perturb", "random", "mix").index(d["start_pool"]))
            dlg.corr_spin.setValue(d["noise_correlation_us"]); dlg.rounds_spin.setValue(d["screen_rounds"])
            dlg.screen_spin.setValue(d["screen_nfev"]); dlg.keep_spin.setValue(d["keep_fraction"])
            dlg.cross_spin.setValue(d["crossover_rate"]); dlg.kicks_spin.setValue(d["n_kicks"])
            dlg.kick_size_spin.setValue(d["kick_size"]); dlg.temp_spin.setValue(d["kick_temperature"])
            dlg.coarse_edit.setText(", ".join(str(m) for m in d["coarse_levels"]))
        finally:
            self._applying = False
        self._preset_tip(name)

    def _preset_changed(self, name):
        if not self._applying:
            if name == "custom":
                self._preset_tip(name)
            else:
                self.apply_preset(name)

    def _mark_custom(self, *args):
        if not self._applying and self.preset_box.currentText() != "custom":
            self._applying = True
            self.preset_box.setCurrentText("custom")
            self._applying = False
            self._preset_tip("custom")

    def _preset_tip(self, name):
        self.preset_box.setToolTip(f"{name}: {PRESET_TIPS[name]}\n(exploitation ↔ exploration; details under Advanced…)")

    def _open_advanced(self):
        self.explore_dialog.show(); self.explore_dialog.raise_(); self.explore_dialog.activateWindow()

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
            multistart_noise=self.spread_spin.value(), reweight_rounds=self.reweight_spin.value(),
            reweight_power=self.reweight_power_spin.value(), limits_late=self.limits_late_check.isChecked(),
            **self._exploration_fields())
        if cfg.controls == "spline" and self.start_pulse is not None and cfg.n_spline > self.start_pulse.N:
            raise ValueError(f"more spline control points ({cfg.n_spline}) than sub-pulses ({self.start_pulse.N})")
        cfg.check()
        return cfg

    def _exploration_fields(self):
        d = self.explore_dialog
        return dict(start_pool=("perturb", "random", "mix")[d.pool_box.currentIndex()],
                    noise_correlation_us=d.corr_spin.value(), screen_rounds=d.rounds_spin.value(),
                    screen_nfev=d.screen_spin.value(), keep_fraction=d.keep_spin.value(),
                    crossover_rate=d.cross_spin.value(), n_kicks=d.kicks_spin.value(), kick_size=d.kick_size_spin.value(),
                    kick_temperature=d.temp_spin.value(), coarse_levels=d.coarse_levels())

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
        if self.preset_box.currentText() in EXPLORATION_PRESETS:
            self.apply_preset(self.preset_box.currentText())
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
                self.smooth_edit, self.amp_rate_spin, self.phase_rate_spin, self.starts_spin, self.spread_spin,
                self.preset_box, self.reweight_spin, self.reweight_power_spin, self.limits_late_check] + self.gen_widgets + self.explore_dialog.widgets

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
        self._cfg = cfg
        n_workers = min(self.max_parallel(), max(int(cfg.n_starts), min(int(cfg.n_kicks), self.max_parallel()), 1))
        try:
            self._search = Search(self.start_pulse, cfg, settings, parallel=n_workers)
        except ValueError as e:
            self._show(str(e), BAD)
            return
        self._limit, self._started = settings.max_nfev, time.time()
        self._clear_convergence()
        self._current = {}
        ctx = mp.get_context("spawn")
        self._tasks, self._queue, self._stop_event = ctx.Queue(), ctx.Queue(), ctx.Event()
        self._workers = [ctx.Process(target=grape_worker, args=(self._tasks, self._queue, self._stop_event), daemon=True)
                         for _ in range(n_workers)]
        for proc in self._workers:
            proc.start()
        self._process = True
        self._dispatch()
        self._stop_requested_at = None
        self._poll = QTimer(self); self._poll.setInterval(100); self._poll.timeout.connect(self._poll_queue)
        self._poll.start()
        self.solve_button.setText("Stop")
        self._set_inputs_enabled(False)
        self._show("starting optimizer …")

    def _dispatch(self):
        for task in self._search.tasks():
            self._tasks.put(task)
        if self._search.done:
            self._end_workers()

    def _end_workers(self):
        for _ in self._workers:
            self._tasks.put(None)

    def stop(self, wait=False):
        if not self.solving:
            return
        self._stop_event.set()
        self._search.stop()
        if self._stop_requested_at is None:
            self._stop_requested_at = time.time()
        self._show("stopping …")
        if wait:
            self._end_workers()
            for proc in self._workers:
                proc.join(5)
                if proc.is_alive():
                    proc.terminate()
            self._finish()

    def _progress_text(self):
        s = self._search
        best = min(self._current.values()) if self._current else float("nan")
        evals = sum(h[0][-1] for h in self._history.values() if h[0])
        if s.phase == "rounds":
            n = int(self._cfg.screen_rounds)
            stage = (f"round {s.round + 1} of {n + 1}" if n and s.round < n else
                     "final round" if n else (f"{len(s.cands)} starts" if len(s.cands) > 1 else ""))
        else:
            stage = f"kicks ({s.kicks_done} of {int(self._cfg.n_kicks)})"
        alive = sum(1 for c in s.cands.values() if c["alive"])
        multi = f"{stage}, {alive} candidates, {s.running} running, " if len(s.cands) > 1 or s.phase != "rounds" else ""
        return (f"optimizing … {multi}{evals} evaluations (limit {self._limit} per candidate), "
                f"{time.time() - self._started:.0f} s, best goal residual {best:.2e}  (click Stop to interrupt)")

    def _set_history(self, cid, n, c):
        self._history[cid] = (n, c)
        self._conv_curve(cid).setData(np.array(n, float), np.maximum(np.array(c, float), 1e-300))

    def _poll_queue(self):
        changed = False
        while True:
            try:
                msg = self._queue.get_nowait()
            except queue_module.Empty:
                break
            except (EOFError, OSError):
                self._search.stop()
                self._failure = "an optimizer process ended unexpectedly"
                break
            kind = msg[0]
            if kind == "progress":
                cid = msg[1]
                n, c = self._history.get(cid, ([], []))
                self._set_history(cid, n + [msg[2]], c + [msg[3]])
                self._current[cid] = msg[3]
                changed = True
            elif kind == "snapshot":
                cid = msg[1]
                self._current[cid] = msg[4]
                if msg[4] <= min(self._current.values()):           # show the candidate that is best right now
                    ps = PulseSequence(np.asarray(msg[2]), np.asarray(msg[3]), self.start_pulse.dt, self.start_pulse.umax)
                    self.result = ps
                    self._show_result(ps)
                    self.pulseChanged.emit(ps)
            elif kind == "done":
                res = msg[1]
                cid = res["cid"]
                info = res.get("info")
                if info is not None:
                    self._current[cid] = info["goal_residual"]
                    hist = list(np.asarray(info.get("history", []), float))
                    if hist:
                        n0 = info["nfev"] - len(hist)
                        n, c = self._history.get(cid, ([], []))
                        keep = [k for k, v in enumerate(n) if v <= n0]
                        self._set_history(cid, [n[k] for k in keep] + list(range(n0 + 1, info["nfev"] + 1)),
                                          [c[k] for k in keep] + hist)
                rnd = (self._search.phase, self._search.round)
                self._search.submit(res)
                if (self._search.phase, self._search.round) != rnd:
                    self._recolor()
                self._dispatch()
                changed = True
        if self._stop_requested_at is not None and time.time() - self._stop_requested_at > 5:
            for proc in self._workers:
                proc.terminate()
            self._search.phase = "done"
        if not self._search.done and not any(p.is_alive() for p in self._workers) and self._queue.empty():
            self._failure = "the optimizer processes ended without a result"
            self._search.phase = "done"
        if self._search.done and (not self._search.running or self._stop_requested_at is not None
                                  or not any(p.is_alive() for p in self._workers)):
            self._complete()
        elif changed:
            self._show(self._progress_text())

    def _recolor(self, best=None):
        for cid, curve in self.conv_curves.items():
            c = self._search.cands.get(cid)
            if c is not None and not c["alive"] and cid != best:
                curve.setPen(DEAD_PEN)
            else:
                curve.setPen(pg.mkPen(START_COLORS[cid % len(START_COLORS)], width=3 if cid == best else 1.5))

    def _complete(self):
        self._finish()
        try:
            cid, pulse, info = self._search.best()
        except RuntimeError as e:
            self._show(getattr(self, "_failure", None) or str(e) or "optimizer terminated (no result)", BAD)
            return
        self.result = pulse
        self._show_result(self.result)
        self._recolor(best=cid)
        self.last_infos = {k: c["info"] for k, c in self._search.cands.items() if c["info"] is not None}
        self._report(info)
        self.pulseChanged.emit(self.result)

    def _finish(self):
        if getattr(self, "_poll", None) is not None:
            self._poll.stop()
            self._poll = None
        if getattr(self, "_workers", None):
            self._end_workers()
        for proc in getattr(self, "_workers", []):
            proc.join(0.5)
            if proc.is_alive():
                proc.terminate()
        self._workers = []
        self._process = None
        self.solve_button.setText("Optimize")
        self._set_inputs_enabled(True)

    def _report(self, info):
        reason = info.get("stop_reason")
        if info["reached"]:
            state = "all goals reached"
        else:
            state = {"user": "stopped – best so far", "limit": "evaluation limit reached – best so far",
                     "time": "time limit reached – best so far", "paused": "stopped – best so far"}.get(
                reason, "converged, goals not fully reached")
        fmt = lambda d: ", ".join(f"{k} = {v:.3g}" for k, v in d.items())
        srch = info.get("search", {})
        n, evals = srch.get("candidates", 1), srch.get("evaluations", info["nfev"])
        if n > 1:
            kind = {"start": "start", "child": "crossover child", "kick": "kick"}.get(srch.get("kind"), "candidate")
            state += f" – best of {n} candidates: {kind} {srch.get('best', 0) + 1} ({srch.get('reached', 0)} reached the goals)"
        log = "\n".join(srch.get("log", [])) if n > 1 else ""
        self._show(f"{state} ({evals} evaluations, {time.time() - self._started:.1f} s)\nresult: {fmt(info['values'])}"
                   f"\nstart: {fmt(info.get('start_summary', {}))}" + (f"\nsearch: {log}" if log else ""),
                   GOOD if info["reached"] else BAD)

    def _show(self, text, color=None):
        self.status.setStyleSheet(f"color: {color};" if color else "")
        self.status.setText(text)
        self._fit_label(self.status)

    @staticmethod
    def _fit_label(label):
        # a word-wrapped label under-reports its height: reserve the wrapped height so that the layout
        # (inside the scroll area) grows instead of squeezing the other inputs
        width = label.width() if label.width() > 50 else 400
        label.setMinimumHeight(label.heightForWidth(width) if label.text() else 0)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._fit_label(self.status)
