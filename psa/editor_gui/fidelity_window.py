"""Fidelity analysis window of the editor (opened with ⊕ on the panel's fidelity plot).

Main plot: fidelity over offset × amplitude (B1 scaling) as a 3D surface or a heat map, with settable
ranges and resolution. Below (scrollable): fidelity vs. offset at a fixed amplitude and fidelity vs.
amplitude at a fixed offset. Reference: the pulse's own gate at 0 kHz / 100 % (as the panel), a UR rotation
or a PP transfer. Not live: "Refresh" takes the editor's current pulse and computes everything in a separate
worker process (psa.design.fidelity), so the window never blocks and does not interfere with a running
optimization; a new Refresh replaces a running one.
"""
from __future__ import annotations

import multiprocessing as mp
import queue as queue_module
import time

import numpy as np
import pyqtgraph as pg
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFormLayout, QGroupBox, QHBoxLayout,
                               QLabel, QPushButton, QScrollArea, QSpinBox, QVBoxLayout, QWidget)

from psa.design.fidelity import fidelity_worker, summary
from psa.editor_gui.widgets import (ACCENT, BAD, GOOD, MUTED_FG, PANEL_BG, TEXT_FG, lock_heights, make_plot,
                                    spawn_single_blas)

QUANTITIES = ("log10(1 − F)", "F")
REFERENCES = ("automatic: the GRAPE gate goal while optimizing there, else the pulse's own gate",
              "pulse's own gate (0 kHz, 100 %)", "gate goal of the GRAPE tab", "UR rotation (angle about an axis)",
              "PP (m0 → target)")
AUTO, OWN, GOAL, UR, PP = range(5)


def _fmt(v):
    return "(" + ", ".join(f"{float(x):g}" for x in v) + ")"


def describe(ref):
    """Short text of a reference dict."""
    if ref.get("kind") == "UR":
        return f"UR {float(ref['angle']):g}° about {_fmt(ref['axis'])}"
    if ref.get("kind") == "PP":
        return f"PP {_fmt(ref['m0'])} → {_fmt(ref['mt'])}"
    if ref.get("kind") == "points":
        n_own = sum(len(p) > 3 and p[3] is not None for p in ref["points"])
        free = {"phase": "; up to a free phase", "linear": "; up to a linear phase (phase + delay)"}.get(
            ref.get("freedom", "none"), "")
        return (f"goals per point ({len(ref['points'])} points, {n_own} with own goal; between the points: the goal "
                f"of the nearest point{free})")
    return "the pulse's own gate (0 kHz, 100 %)"


def _spin(lo, hi, value, decimals=2, suffix=""):
    sp = QDoubleSpinBox(); sp.setRange(lo, hi); sp.setDecimals(decimals); sp.setValue(value); sp.setSuffix(suffix)
    return sp


def _count(value, lo=2, hi=2001, suffix=" points"):
    sp = QSpinBox(); sp.setRange(lo, hi); sp.setValue(value); sp.setSuffix(suffix)
    return sp


def _vector(values):
    row = QHBoxLayout(); spins = []
    for name, v in zip("xyz", values):
        sp = _spin(-1, 1, v, 3); sp.setPrefix(name + " "); sp.setSingleStep(0.1)
        row.addWidget(sp); spins.append(sp)
    return row, spins


def _range_row(lo, hi, n, unit, rng):
    a, b = _spin(*rng, lo, 2, unit), _spin(*rng, hi, 2, unit)
    c = _count(n)
    row = QHBoxLayout(); row.addWidget(a); row.addWidget(QLabel("to")); row.addWidget(b); row.addWidget(c)
    return row, (a, b, c)


class FidelityWindow(QDialog):
    def __init__(self, pulse_provider, parent=None, goal_provider=None):
        super().__init__(parent)
        self.setWindowTitle("Fidelity analysis")
        self.setWindowFlag(Qt.WindowType.Window, True)
        self.resize(900, 900)
        self.pulse_provider = pulse_provider          # callable() -> PulseSequence or None
        self.goal_provider = goal_provider            # callable() -> (GRAPE gate goal or None, active)
        self.result = None                            # last result (+ the settings it was computed with)
        self._computed_pulse = None
        self._pending = None                          # id of the running job
        self._job_id = 0
        self._proc = None
        self._started = None

        outer = QVBoxLayout(self)
        top = QHBoxLayout()
        self.refresh_button = QPushButton("Refresh"); self.refresh_button.clicked.connect(self.refresh)
        self.refresh_button.setToolTip("compute all enabled plots for the editor's current pulse (in a separate "
                                       "process; a running computation is replaced)")
        self.status = QLabel("not computed yet: click Refresh"); self.status.setWordWrap(True)
        self.status.setStyleSheet(f"color: {MUTED_FG};")
        top.addWidget(self.refresh_button); top.addWidget(self.status, 1)
        outer.addLayout(top)
        sf = QFormLayout()
        self.quantity_box = QComboBox(); self.quantity_box.addItems(QUANTITIES)
        self.quantity_box.currentIndexChanged.connect(self.redraw)
        sf.addRow("show:", self.quantity_box)
        self.ref_box = QComboBox(); self.ref_box.addItems(REFERENCES)
        self.ref_box.setToolTip("what the pulse should do. automatic: the gate goal of the GRAPE tab while it optimizes "
                                "or is shown, else the pulse's own gate at 0 kHz / 100 % (as the panel's plot); or "
                                "fixed: own gate, the GRAPE gate goal, a UR rotation or a PP transfer (set below). "
                                "F = (3 + tr(Rt^T R))/6 (UR, own) or (1 + R m0 · target)/2 (PP); taken at Refresh")
        sf.addRow("reference:", self.ref_box)
        self.ref_info = QLabel(""); self.ref_info.setStyleSheet(f"color: {MUTED_FG};")
        sf.addRow("", self.ref_info)
        self.ur_widget = QWidget(); uf = QFormLayout(self.ur_widget); uf.setContentsMargins(0, 0, 0, 0)
        self.angle_spin = _spin(-720, 720, 180.0, 3, " °")
        row, self.axis_spins = _vector((1.0, 0.0, 0.0))
        uf.addRow("angle:", self.angle_spin); uf.addRow("axis:", row)
        self.pp_widget = QWidget(); pf = QFormLayout(self.pp_widget); pf.setContentsMargins(0, 0, 0, 0)
        row, self.m0_spins = _vector((0.0, 0.0, 1.0)); pf.addRow("from m0:", row)
        row, self.mt_spins = _vector((1.0, 0.0, 0.0)); pf.addRow("to:", row)
        sf.addRow(self.ur_widget); sf.addRow(self.pp_widget)
        self._ref_edited = False
        for sp in [self.angle_spin] + self.axis_spins + self.m0_spins + self.mt_spins:
            sp.valueChanged.connect(self._mark_ref_edited)
        self.ref_box.currentIndexChanged.connect(self._ref_changed)
        outer.addLayout(sf)

        content = QWidget(); cl = QVBoxLayout(content)
        # ---- main map
        mbox = QGroupBox("Fidelity vs. offset and amplitude")
        ml = QVBoxLayout(mbox)
        mf = QFormLayout()
        self.map_check = QCheckBox("compute on Refresh"); self.map_check.setChecked(True)
        mf.addRow(self.map_check)
        row, self.map_off = _range_row(-10.0, 10.0, 101, " kHz", (-1e4, 1e4)); mf.addRow("offset:", row)
        row, self.map_amp = _range_row(80.0, 120.0, 41, " %", (0.0, 1e4)); mf.addRow("amplitude:", row)
        self.view_box = QComboBox(); self.view_box.addItems(["surface (3D)", "heat map"])
        self.view_box.currentIndexChanged.connect(self.redraw)
        mf.addRow("view:", self.view_box)
        ml.addLayout(mf)
        self.figure = Figure(figsize=(7, 5)); self.figure.patch.set_facecolor(PANEL_BG)
        self.canvas = FigureCanvas(self.figure); self.canvas.setMinimumHeight(480)
        ml.addWidget(self.canvas)
        self.map_info = QLabel(""); self.map_info.setWordWrap(True)
        ml.addWidget(self.map_info)
        cl.addWidget(mbox)
        # ---- cut: vs offset at a fixed amplitude
        obox = QGroupBox("Fidelity vs. offset at a fixed amplitude")
        ol = QVBoxLayout(obox); of = QFormLayout()
        self.cut_off_check = QCheckBox("compute on Refresh"); self.cut_off_check.setChecked(True)
        of.addRow(self.cut_off_check)
        self.cut_off_amp = _spin(0.0, 1e4, 100.0, 2, " %"); of.addRow("amplitude:", self.cut_off_amp)
        row, self.cut_off_range = _range_row(-15.0, 15.0, 601, " kHz", (-1e4, 1e4)); of.addRow("offset:", row)
        ol.addLayout(of)
        self.cut_off_plot = make_plot("log10(1 − F)", "offset [kHz]"); self.cut_off_plot.setMinimumHeight(260)
        self.cut_off_curve = self.cut_off_plot.plot(pen=pg.mkPen(ACCENT, width=2), connect="finite")
        ol.addWidget(self.cut_off_plot)
        self.cut_off_info = QLabel(""); ol.addWidget(self.cut_off_info)
        cl.addWidget(obox)
        # ---- cut: vs amplitude at a fixed offset
        abox = QGroupBox("Fidelity vs. amplitude at a fixed offset")
        al = QVBoxLayout(abox); af = QFormLayout()
        self.cut_amp_check = QCheckBox("compute on Refresh"); self.cut_amp_check.setChecked(True)
        af.addRow(self.cut_amp_check)
        self.cut_amp_off = _spin(-1e4, 1e4, 0.0, 3, " kHz"); af.addRow("offset:", self.cut_amp_off)
        row, self.cut_amp_range = _range_row(50.0, 150.0, 201, " %", (0.0, 1e4)); af.addRow("amplitude:", row)
        al.addLayout(af)
        self.cut_amp_plot = make_plot("log10(1 − F)", "amplitude [%]"); self.cut_amp_plot.setMinimumHeight(260)
        self.cut_amp_curve = self.cut_amp_plot.plot(pen=pg.mkPen(ACCENT, width=2), connect="finite")
        al.addWidget(self.cut_amp_plot)
        self.cut_amp_info = QLabel(""); al.addWidget(self.cut_amp_info)
        cl.addWidget(abox)
        cl.addStretch()
        scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(content)
        outer.addWidget(scroll, 1)

        self._poll = QTimer(self); self._poll.setInterval(100); self._poll.timeout.connect(self._poll_results)
        self._ref_changed()
        lock_heights(self)
        self.redraw()

    # ------------------------------------------------------------ inputs ---
    def _mark_ref_edited(self, *args):
        if not getattr(self, "_filling", False):
            self._ref_edited = True

    def _ref_changed(self, *args):
        i = self.ref_box.currentIndex()
        if i in (UR, PP) and not getattr(self, "_ref_edited", False):
            goal = self.goal_provider()[0] if self.goal_provider else None
            if goal is not None and goal.get("kind") != "points":                      # start from the GRAPE gate goal
                self.set_reference_inputs(goal)
        self.ur_widget.setVisible(i == UR)
        self.pp_widget.setVisible(i == PP)
        self.update_ref_info()

    def update_ref_info(self):
        """Show which reference the current choice means right now."""
        if self.ref_box.currentIndex() in (AUTO, GOAL):
            try:
                self.ref_info.setText("now: " + self.reference()[1])
            except ValueError as e:
                self.ref_info.setText(str(e))
        else:
            self.ref_info.setText("")

    def showEvent(self, event):
        super().showEvent(event)
        self.update_ref_info()

    def set_reference_inputs(self, ref):
        """Fill the UR / PP inputs from a reference dict."""
        self._filling = True
        try:
            self._fill(ref)
        finally:
            self._filling = False

    def _fill(self, ref):
        if ref.get("kind") == "UR":
            self.angle_spin.setValue(float(ref["angle"]))
            for sp, v in zip(self.axis_spins, ref["axis"]):
                sp.setValue(float(v))
        elif ref.get("kind") == "PP":
            for sp, v in zip(self.m0_spins + self.mt_spins, list(ref["m0"]) + list(ref["mt"])):
                sp.setValue(float(v))

    def reference(self):
        """(reference dict, description) of the current choice (ValueError if the GRAPE goal is wanted but unset)."""
        i = self.ref_box.currentIndex()
        if i in (AUTO, GOAL):
            goal, active = self.goal_provider() if self.goal_provider else (None, False)
            if i == GOAL and goal is None:
                raise ValueError("the GRAPE tab has no gate goal (UR or PP)")
            if goal is not None and (i == GOAL or active):
                return goal, "GRAPE gate goal: " + describe(goal)
            return {"kind": "own"}, "the pulse's own gate (0 kHz, 100 %)"
        if i == UR:
            ref = {"kind": "UR", "axis": [s.value() for s in self.axis_spins], "angle": self.angle_spin.value()}
        elif i == PP:
            ref = {"kind": "PP", "m0": [s.value() for s in self.m0_spins], "mt": [s.value() for s in self.mt_spins]}
        else:
            return {"kind": "own"}, "the pulse's own gate (0 kHz, 100 %)"
        return ref, describe(ref)

    @staticmethod
    def _values(spins, scale=1.0):
        a, b, n = spins
        if not b.value() > a.value():
            raise ValueError("each range needs 'to' above 'from'")
        return np.linspace(a.value(), b.value(), n.value()) * scale

    def job(self, ps):
        """The job for a pulse from the inputs (ValueError if invalid)."""
        ref, label = self.reference()
        for key in ("axis", "m0", "mt"):
            if ref["kind"] != "points" and key in ref and not any(abs(v) > 0 for v in ref[key]):
                raise ValueError("the reference vectors must not be zero")
        job = {"pulse": [ps.amplitude.tolist(), ps.phase.tolist(), ps.dt, ps.umax], "reference": ref,
               "reference_label": label,
               "map": None, "cut_offset": None, "cut_amp": None}
        if self.map_check.isChecked():
            scales = self._values(self.map_amp, 0.01)
            if scales[0] <= 0:
                raise ValueError("amplitudes must be positive")
            job["map"] = {"offsets_hz": self._values(self.map_off, 1e3).tolist(), "scales": scales.tolist()}
        if self.cut_off_check.isChecked():
            if self.cut_off_amp.value() <= 0:
                raise ValueError("the fixed amplitude must be positive")
            job["cut_offset"] = {"offsets_hz": self._values(self.cut_off_range, 1e3).tolist(),
                                 "scale": self.cut_off_amp.value() / 100}
        if self.cut_amp_check.isChecked():
            scales = self._values(self.cut_amp_range, 0.01)
            if scales[0] <= 0:
                raise ValueError("amplitudes must be positive")
            job["cut_amp"] = {"scales": scales.tolist(), "offset_hz": self.cut_amp_off.value() * 1e3}
        if not any(job[k] for k in ("map", "cut_offset", "cut_amp")):
            raise ValueError("nothing to compute: enable at least one plot")
        return job

    # ------------------------------------------------------------- worker ---
    @property
    def busy(self):
        return self._pending is not None

    def _ensure_worker(self):
        if self._proc is None or not self._proc.is_alive():
            ctx = mp.get_context("spawn")
            self._jobs, self._results = ctx.Queue(), ctx.Queue()
            self._proc = spawn_single_blas(ctx, fidelity_worker, (self._jobs, self._results))[0]

    def _kill_worker(self):
        if self._proc is not None:
            if self._proc.is_alive():
                self._proc.terminate()
            self._proc.join(1)
            self._proc = None

    def refresh(self):
        self.update_ref_info()
        ps = self.pulse_provider() if self.pulse_provider else None
        if ps is None or ps.N == 0:
            self._show("no pulse in the editor", BAD)
            return
        try:
            job = self.job(ps)
        except ValueError as e:
            self._show(str(e), BAD)
            return
        if self.busy:                                 # replace the running computation
            self._kill_worker()
        self._job_id += 1
        job["id"] = self._job_id
        self._job = job
        self._job_pulse = ps
        self._ensure_worker()
        self._jobs.put(job)
        self._pending, self._started = self._job_id, time.time()
        self._poll.start()
        self._show(f"computing for {ps.N} sub-pulses, {ps.duration * 1e6:.1f} µs …")

    def _poll_results(self):
        while True:
            try:
                res = self._results.get_nowait()
            except queue_module.Empty:
                break
            except (EOFError, OSError):
                self._pending = None
                self._show("the computation process ended unexpectedly", BAD)
                break
            if res.get("id") != self._pending:
                continue                              # a replaced job
            self._pending = None
            if "error" in res:
                self._show(f"computation failed: {res['error']}", BAD)
            else:
                self.result = dict(res, job=self._job)
                self._computed_pulse = self._job_pulse
                self._show(f"computed in {res['elapsed_s']:.1f} s for {self._job_pulse.N} sub-pulses, "
                           f"{self._job_pulse.duration * 1e6:.1f} µs ({time.strftime('%H:%M:%S')}); reference: "
                           f"{self._job['reference_label']}", GOOD)
                self.redraw()
        if self._pending is None:
            self._poll.stop()
        elif self._proc is not None and not self._proc.is_alive():
            self._pending = None
            self._poll.stop()
            self._show("the computation process ended unexpectedly", BAD)

    def pulse_changed(self, ps):
        """The editor's pulse changed: mark the plots as outdated."""
        if self._computed_pulse is not None and ps is not self._computed_pulse and not self.busy:
            self.status.setText(self.status.text().split("  —")[0] + "  — the pulse has changed since: Refresh")
            self.status.setStyleSheet(f"color: {MUTED_FG};")

    def _show(self, text, color=None):
        self.status.setText(text)
        self.status.setStyleSheet(f"color: {color or MUTED_FG};")

    # -------------------------------------------------------------- plots ---
    def _y(self, inf):
        inf = np.maximum(np.asarray(inf, float), 0.0)
        return np.log10(np.maximum(inf, 1e-16)) if self.quantity_box.currentIndex() == 0 else 1 - inf

    @staticmethod
    def _stats(inf):
        s = summary(inf)
        return (f"worst 1 − F {s['worst']:.2e}, mean {s['mean']:.2e}; below 1e-3: {100 * s['below 0.001']:.0f} %, "
                f"1e-4: {100 * s['below 0.0001']:.0f} %, 1e-5: {100 * s['below 1e-05']:.0f} %")

    def redraw(self, *args):
        label = QUANTITIES[self.quantity_box.currentIndex()]
        r = self.result
        self._draw_map(r, label)
        for plot, curve, info, key, xs in (
                (self.cut_off_plot, self.cut_off_curve, self.cut_off_info, "cut_offset",
                 lambda j: np.asarray(j["cut_offset"]["offsets_hz"]) / 1e3),
                (self.cut_amp_plot, self.cut_amp_curve, self.cut_amp_info, "cut_amp",
                 lambda j: np.asarray(j["cut_amp"]["scales"]) * 100)):
            plot.setLabel("left", label)
            if r is not None and r.get(key) is not None:
                curve.setData(xs(r["job"]), self._y(r[key]))
                info.setText(self._stats(r[key]))
            else:
                curve.setData([], []); info.setText("")
        if r is not None and r["job"].get("cut_offset"):
            self.cut_off_plot.setTitle(f"at {r['job']['cut_offset']['scale'] * 100:g} %", color=MUTED_FG, size="9pt")
        if r is not None and r["job"].get("cut_amp"):
            self.cut_amp_plot.setTitle(f"at {r['job']['cut_amp']['offset_hz'] / 1e3:g} kHz", color=MUTED_FG, size="9pt")

    def _draw_map(self, r, label):
        fig = self.figure
        fig.clear()
        if r is None or r.get("map") is None:
            ax = fig.add_subplot(111); ax.set_facecolor(PANEL_BG); ax.set_axis_off()
            ax.text(0.5, 0.5, "click Refresh", color=MUTED_FG, ha="center", va="center", transform=ax.transAxes)
            self.map_info.setText("")
            self.canvas.draw_idle()
            return
        job = r["job"]["map"]
        o = np.asarray(job["offsets_hz"]) / 1e3
        a = np.asarray(job["scales"]) * 100
        Z = self._y(r["map"])                        # (n_amp, n_off)
        if self.view_box.currentIndex() == 0:
            ax = fig.add_subplot(111, projection="3d")
            X, Y = np.meshgrid(o, a)
            surf = ax.plot_surface(X, Y, Z, cmap="jet", linewidth=0, antialiased=False, rcount=min(len(a), 120),
                                   ccount=min(len(o), 160))
            ax.set_zlabel(label, color=TEXT_FG, fontsize=8)
            for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
                axis.set_pane_color((0.07, 0.09, 0.15, 1.0))
            cb = fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.1)
        else:
            ax = fig.add_subplot(111)
            im = ax.imshow(Z, origin="lower", aspect="auto", cmap="jet",
                           extent=(o[0] - (o[1] - o[0]) / 2, o[-1] + (o[1] - o[0]) / 2,
                                   a[0] - (a[1] - a[0]) / 2, a[-1] + (a[1] - a[0]) / 2))
            cut = r["job"]
            if cut.get("cut_offset"):
                ax.axhline(cut["cut_offset"]["scale"] * 100, color="white", lw=0.8, ls="--")
            if cut.get("cut_amp"):
                ax.axvline(cut["cut_amp"]["offset_hz"] / 1e3, color="white", lw=0.8, ls="--")
            cb = fig.colorbar(im, ax=ax)
            cb.set_label(label, color=TEXT_FG)
        ax.set_facecolor(PANEL_BG)
        ax.set_xlabel("offset [kHz]", color=TEXT_FG, fontsize=9)
        ax.set_ylabel("amplitude [%]", color=TEXT_FG, fontsize=9)
        ax.tick_params(colors=MUTED_FG, labelsize=8)
        cb.ax.tick_params(colors=MUTED_FG, labelsize=8)
        if self.view_box.currentIndex() == 0:
            fig.subplots_adjust(left=0.0, right=0.95, bottom=0.03, top=0.99)
        else:
            fig.subplots_adjust(left=0.1, right=0.97, bottom=0.12, top=0.96)
        self.map_info.setText(self._stats(r["map"]))
        self.canvas.draw_idle()

    def closeEvent(self, event):
        self._poll.stop()
        self._pending = None
        if self._proc is not None:
            try:
                self._jobs.put(None)
            except (OSError, ValueError):
                pass
            self._kill_worker()
        super().closeEvent(event)
