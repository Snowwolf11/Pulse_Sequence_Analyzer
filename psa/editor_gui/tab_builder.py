"""Curve builder tab: compose the error curve from segments (line, arc, helix, Bezier) with twists,
see it in 3D and in three projections, drag Bezier control points, and solve for closure /
vanishing areas / third order."""
from __future__ import annotations

import json
import math

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (QCheckBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGridLayout, QGroupBox,
                               QHBoxLayout, QLabel, QListWidget, QPushButton, QSplitter, QVBoxLayout, QWidget)

from psa.design.builder import PARAMS, UNITS, CurveDesign, Segment
from psa.editor_gui.solver_box import SolverBox
from psa.editor_gui.widgets import ACCENT, ACCENT_2, BAD, GOOD, MUTED_FG, KeyframeItem, Mpl3D, lock_heights, make_plot

PROJECTIONS = (("xy", 0, 1), ("xz", 0, 2), ("yz", 1, 2))
LABELS = {"length": "length", "radius": "radius", "pitch": "pitch", "angle": "angle", "twist": "twist (phase jump)",
          "p1": "P1 (along tangent)", "p2x": "P2 x'", "p2y": "P2 y'", "p2z": "P2 t", "p3x": "P3 x'", "p3y": "P3 y'",
          "p3z": "P3 t"}


def describe(i, seg):
    p = seg.params
    text = {"line": lambda: f"line {p['length']:.4g} µs",
            "arc": lambda: f"arc R={p['radius']:.4g} µs, {p['angle']:.4g}°",
            "helix": lambda: f"helix ρ={p['radius']:.4g} µs, h={p['pitch']:.4g} µs, {p['angle']:.4g}°",
            "bezier": lambda: "Bézier"}[seg.kind]()
    twist = f", twist {seg.twist:.4g}°" if seg.twist else ""
    free = f"  [free: {', '.join(seg.free)}]" if seg.free else ""
    return f"{i + 1}. {text}{twist}{free}"


class BuilderTab(QWidget):
    pulseChanged = Signal(object)
    editInWaveform = Signal(object)

    def __init__(self, dt=0.5e-6, umax=1e4, parent=None):
        super().__init__(parent)
        rmin = 1e6 / (2 * math.pi * umax)
        self.design = CurveDesign([Segment("arc", {"radius": 1.2 * rmin, "angle": 180.0}),
                                   Segment("arc", {"radius": 1.2 * rmin, "angle": 180.0}, twist=90.0)], umax, dt)
        self.pulse = None
        self.editors = {}
        outer = QHBoxLayout(self)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter)

        # ---------------------------------------------------------- left side
        left = QWidget(); ll = QVBoxLayout(left)
        self.left_widget = left
        seg_group = QGroupBox("Segments")
        seg_group.setToolTip("each segment starts where the previous one ends, in the same direction")
        sl = QVBoxLayout(seg_group)
        self.seg_list = QListWidget(); self.seg_list.setMinimumHeight(140)
        sl.addWidget(self.seg_list)
        grid = QGridLayout()
        buttons = [("+ Line", lambda: self.add_segment("line")), ("+ Arc", lambda: self.add_segment("arc")),
                   ("+ Helix", lambda: self.add_segment("helix")), ("+ Bézier", lambda: self.add_segment("bezier")),
                   ("Duplicate", self.duplicate_segment), ("Remove", self.remove_segment),
                   ("▲ up", lambda: self.move_segment(-1)), ("▼ down", lambda: self.move_segment(1))]
        for k, (text, slot) in enumerate(buttons):
            b = QPushButton(text); b.clicked.connect(slot)
            grid.addWidget(b, k // 3, k % 3)
        sl.addLayout(grid)
        prow = QHBoxLayout()
        prow.addWidget(QLabel("start phase:"))
        self.phase0_spin = QDoubleSpinBox(); self.phase0_spin.setRange(-1e5, 1e5); self.phase0_spin.setSuffix(" °")
        self.phase0_spin.valueChanged.connect(self._phase0_changed)
        prow.addWidget(self.phase0_spin); prow.addStretch()
        sl.addLayout(prow)
        ll.addWidget(seg_group)

        self.edit_group = QGroupBox("Selected segment")
        self.edit_form = QFormLayout(self.edit_group)
        ll.addWidget(self.edit_group)
        self.rmin_label = QLabel(""); self.rmin_label.setStyleSheet(f"color: {MUTED_FG};"); self.rmin_label.setWordWrap(True)
        ll.addWidget(self.rmin_label)

        self.solver_box = SolverBox("Solve", show_variables=False,
                                    variables_hint="free variables: the segment parameters ticked 'free'")
        self.solver_box.job_provider = lambda goals, keys: {"kind": "builder", "design": self.design.to_dict(), "free": None}
        self.solver_box.solved.connect(self._apply_solution)
        self.solver_box.running.connect(lambda on: self._set_editing_enabled(not on))
        # names kept for scripts / tests
        self.solve_button, self.order_checks, self.solve_label = (self.solver_box.button, self.solver_box.order_checks,
                                                                  self.solver_box.status)
        ll.addWidget(self.solver_box)

        self.info = QLabel(""); self.info.setWordWrap(True)
        ll.addWidget(self.info)
        self.warn = QLabel(""); self.warn.setWordWrap(True); self.warn.setStyleSheet(f"color: {BAD};")
        ll.addWidget(self.warn)
        brow = QHBoxLayout()
        b_edit = QPushButton("Edit in Waveform tab"); b_edit.clicked.connect(lambda: self.pulse is not None and self.editInWaveform.emit(self.pulse))
        b_open = QPushButton("Open design…"); b_open.clicked.connect(self._open_design)
        b_save = QPushButton("Save design…"); b_save.clicked.connect(self._save_design)
        for b in (b_edit, b_open, b_save):
            brow.addWidget(b)
        ll.addLayout(brow)
        ll.addStretch()
        splitter.addWidget(left)

        # --------------------------------------------------------- right side
        right = QWidget(); rl = QGridLayout(right)
        self.view3d = Mpl3D(); self.view3d.setMinimumHeight(300); self.view3d.setMinimumWidth(160)
        rl.addWidget(self.view3d, 0, 0)
        self.proj_plots, self.proj_curves, self.proj_sel, self.proj_marks, self.proj_handles = {}, {}, {}, {}, {}
        positions = {"xy": (0, 1), "xz": (1, 0), "yz": (1, 1)}
        for name, i, j in PROJECTIONS:
            p = make_plot(f"{'xyz'[j]} [µs]", f"{'xyz'[i]} [µs]")
            p.setAspectLocked(True)
            p.setMinimumHeight(260); p.setMinimumWidth(160)
            p.setTitle(f"{name} projection")
            self.proj_curves[name] = p.plot(pen=pg.mkPen(ACCENT, width=1.5))
            self.proj_sel[name] = p.plot(pen=pg.mkPen(ACCENT_2, width=3))
            self.proj_marks[name] = p.plot(pen=None, symbol="o", symbolSize=6, symbolBrush=pg.mkBrush(MUTED_FG))
            handles = KeyframeItem()
            p.addItem(handles)
            handles.sigDragged.connect(lambda idx, x, y, n=name: self._handle_dragged(n, idx, x, y))
            self.proj_handles[name] = handles
            self.proj_plots[name] = p
            rl.addWidget(p, *positions[name])
        hint = QLabel("orange: selected segment · Bézier: drag its control points P1–P3 in the projections")
        hint.setStyleSheet(f"color: {MUTED_FG};"); hint.setWordWrap(True)
        rl.addWidget(hint, 2, 0, 1, 2)
        splitter.addWidget(right)
        splitter.setStretchFactor(1, 2)

        self._timer = QTimer(self); self._timer.setSingleShot(True); self._timer.setInterval(150)
        self._timer.timeout.connect(self.recompute)
        self.seg_list.currentRowChanged.connect(self._selection_changed)
        self._refresh_list(select=0)

    # ------------------------------------------------------------ segments ---
    def selected(self):
        return self.seg_list.currentRow()

    def _refresh_list(self, select=None):
        row = self.selected() if select is None else select
        self.seg_list.blockSignals(True)
        self.seg_list.clear()
        for i, s in enumerate(self.design.segments):
            self.seg_list.addItem(describe(i, s))
        self.seg_list.blockSignals(False)
        if self.design.segments:
            self.seg_list.setCurrentRow(min(max(row, 0), len(self.design.segments) - 1))
        self._build_editor()
        self.schedule()

    def add_segment(self, kind):
        rmin = self.design.r_min_us
        defaults = {"arc": {"radius": 1.2 * rmin}, "helix": {"radius": 1.2 * rmin, "pitch": 2 * rmin}}.get(kind, {})
        i = self.selected() + 1 if self.design.segments else 0
        self.design.segments.insert(i, Segment(kind, defaults))
        self._refresh_list(select=i)

    def duplicate_segment(self):
        i = self.selected()
        if 0 <= i < len(self.design.segments):
            self.design.segments.insert(i + 1, Segment.from_dict(self.design.segments[i].to_dict()))
            self._refresh_list(select=i + 1)

    def remove_segment(self):
        i = self.selected()
        if 0 <= i < len(self.design.segments):
            del self.design.segments[i]
            self._refresh_list(select=i - 1)

    def move_segment(self, step):
        i = self.selected()
        j = i + step
        if 0 <= i < len(self.design.segments) and 0 <= j < len(self.design.segments):
            s = self.design.segments
            s[i], s[j] = s[j], s[i]
            self._refresh_list(select=j)

    def _selection_changed(self, row):
        self._build_editor()
        self._draw_selection()

    def _build_editor(self):
        while self.edit_form.rowCount():
            self.edit_form.removeRow(0)
        self.editors = {}
        i = self.selected()
        if not (0 <= i < len(self.design.segments)):
            self.edit_group.setTitle("Selected segment: none")
            return
        seg = self.design.segments[i]
        self.edit_group.setTitle(f"Selected segment {i + 1}: {seg.kind}")
        for name in [n for n, _ in PARAMS[seg.kind]] + ["twist"]:
            spin = QDoubleSpinBox(); spin.setRange(-1e6, 1e6); spin.setDecimals(4)
            spin.setSuffix(" " + UNITS[name]); spin.setValue(seg.get(name))
            spin.setKeyboardTracking(False)
            spin.valueChanged.connect(lambda v, n=name: self._param_edited(n, v))
            free = QCheckBox("free"); free.setChecked(name in seg.free)
            free.setToolTip("the solver may change this parameter")
            free.toggled.connect(lambda on, n=name: self._free_toggled(n, on))
            row = QHBoxLayout(); row.addWidget(spin, stretch=1); row.addWidget(free)
            w = QWidget(); w.setLayout(row); row.setContentsMargins(0, 0, 0, 0)
            self.edit_form.addRow(LABELS[name] + ":", w)
            self.editors[name] = (spin, free)
        lock_heights(self.edit_group)
        self.rmin_label.setText(f"Tightest radius at full amplitude: R_min = 1/(2π·ν_max) = {self.design.r_min_us:.4g} µs. "
                                "Arc amplitude = R_min/R; helix amplitude = R_min·ρ/(ρ²+c²), c = h/2π. "
                                "Twist = phase jump = rotation of the next segment about the tangent.")

    def _param_edited(self, name, value):
        i = self.selected()
        if 0 <= i < len(self.design.segments):
            self.design.segments[i].set(name, value)
            self.seg_list.item(i).setText(describe(i, self.design.segments[i]))
            self.schedule()

    def _free_toggled(self, name, on):
        i = self.selected()
        if 0 <= i < len(self.design.segments):
            free = self.design.segments[i].free
            if on and name not in free:
                free.append(name)
            elif not on and name in free:
                free.remove(name)
            self.seg_list.item(i).setText(describe(i, self.design.segments[i]))

    def _phase0_changed(self, v):
        self.design.phase0 = v
        self.schedule()

    # -------------------------------------------------------------- compute ---
    def set_metadata(self, dt, umax):
        self.design.dt, self.design.umax = dt, umax
        self._build_editor()
        self.schedule()

    def schedule(self, *args):
        self._timer.start()

    def recompute(self):
        d = self.design
        self.pieces, warnings = d.pieces()
        self.frames = list(d.segment_frames)
        self.points, starts = d.curve(pieces=self.pieces)
        self.starts = dict(starts)
        self.warn.setText("\n".join(warnings))
        if not self.pieces:
            self.info.setText("empty design")
            self.pulse = None
            self._draw()
            self.pulseChanged.emit(None)
            return None
        rb = d.robustness(self.pieces)
        T = d.duration(self.pieces)
        self.pulse = d.to_pulse()
        self.info.setText(
            f"continuous design: {T * 1e6:.3f} µs; pulse: {self.pulse.N} sub-pulses ({self.pulse.duration * 1e6:.3f} µs).\n"
            f"design robustness |A1|/T = {rb['A1_norm']:.2e}, |A2|/T² = {rb['A2_norm']:.2e}, |A3|/T³ = {rb['A3_norm']:.2e} "
            "(continuous; the panel shows the discretized pulse)")
        self._draw()
        self.pulseChanged.emit(self.pulse)
        return self.pulse

    def current_pulse(self):
        return self.pulse if self.pulse is not None else self.recompute()

    @property
    def solving(self):
        return self.solver_box.solving

    def solve(self):
        self.solver_box.toggle()

    def stop_solver(self, wait=False):
        self.solver_box.stop(wait)

    def _set_editing_enabled(self, on):
        from PySide6.QtWidgets import QAbstractButton, QAbstractSpinBox, QAbstractItemView
        inside_box = set(self.solver_box.findChildren(QWidget))
        for cls in (QAbstractButton, QAbstractSpinBox, QAbstractItemView):
            for w in self.left_widget.findChildren(cls):
                if w not in inside_box:
                    w.setEnabled(on)
        for h in self.proj_handles.values():
            h.setEnabled(on)

    def _apply_solution(self, solution, info):
        self.design = CurveDesign.from_dict(solution)
        self.solver_box.status.setText(self.solver_box.status.text() + "\nfree values: " + ", ".join(
            f"seg {i + 1} {n} = {self.design.segments[i].get(n):.5g}" for i, n in self.design.free_parameters()))
        self._refresh_list()
        self._timer.stop()
        self.recompute()

    # -------------------------------------------------------------- drawing ---
    def _segment_range(self, i):
        if i not in self.starts:
            return None
        a = self.starts[i]
        later = [v for k, v in self.starts.items() if k > i]
        b = min(later) if later else len(self.points) - 1
        return a, b

    def _draw(self):
        pts = self.points * 1e6
        self.view3d.plot_curve(self.points, "error curve (design)")
        ax = self.view3d.ax
        rng = self._segment_range(self.selected())
        if rng:
            a, b = rng
            sel = pts[a:b + 1]
            ax.plot(sel[:, 0], sel[:, 1], sel[:, 2], color=ACCENT_2, lw=2.5)
            self.view3d.canvas.draw_idle()
        marks = pts[[v for v in self.starts.values()]] if self.starts else np.zeros((0, 3))
        for name, i, j in PROJECTIONS:
            self.proj_curves[name].setData(pts[:, i], pts[:, j])
            self.proj_marks[name].setData(marks[:, i], marks[:, j])
        self._draw_selection()

    def _draw_selection(self):
        if not hasattr(self, "points"):
            return
        pts = self.points * 1e6
        rng = self._segment_range(self.selected())
        handles = self._control_points()
        for name, i, j in PROJECTIONS:
            if rng:
                a, b = rng
                self.proj_sel[name].setData(pts[a:b + 1, i], pts[a:b + 1, j])
            else:
                self.proj_sel[name].setData([], [])
            if handles is not None:
                self.proj_handles[name].set_points(handles[:, i] * 1e6, handles[:, j] * 1e6)
            else:
                self.proj_handles[name].set_points(np.zeros(0), np.zeros(0))

    def _control_points(self):
        """Lab positions (s) of P1..P3 of the selected Bezier segment, else None."""
        i = self.selected()
        if not (0 <= i < len(self.design.segments)) or self.design.segments[i].kind != "bezier":
            return None
        if i >= len(self.frames):
            return None
        p = self.design.segments[i].params
        local = np.array([[0, 0, p["p1"]], [p["p2x"], p["p2y"], p["p2z"]], [p["p3x"], p["p3y"], p["p3z"]]]) * 1e-6
        R, phase = self.frames[i]
        axes = CurveDesign.local_axes(R, phase)
        start = self.points[self.starts[i]] if i in self.starts else (self.points[-1] if len(self.points) else np.zeros(3))
        self._bezier_frame = (start, axes)
        return local @ axes + start

    def _handle_dragged(self, proj, index, x, y):
        i = self.selected()
        handles = self._control_points()
        if handles is None:
            return
        start, axes = self._bezier_frame
        _, a, b = next(p for p in PROJECTIONS if p[0] == proj)
        lab = handles[index].copy()
        lab[a], lab[b] = x * 1e-6, y * 1e-6
        local = axes @ (lab - start) * 1e6                       # rows of axes are orthonormal
        seg = self.design.segments[i]
        if index == 0:
            seg.set("p1", max(0.0, local[2]))                   # P1 stays on the tangent
        else:
            k = "p2" if index == 1 else "p3"
            for c, v in zip("xyz", local):
                seg.set(k + c, v)
        self._build_editor()
        self._timer.stop()
        self.recompute()

    # -------------------------------------------------------------- designs ---
    def load_design(self, d):
        self.design = CurveDesign.from_dict(d)
        self.phase0_spin.blockSignals(True); self.phase0_spin.setValue(self.design.phase0); self.phase0_spin.blockSignals(False)
        self._refresh_list(select=0)
        self._timer.stop()
        self.recompute()

    def _open_design(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open curve design", "", "Design (*.json)")
        if path:
            with open(path) as fh:
                self.load_design(json.load(fh))

    def _save_design(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save curve design", "curve_design.json", "Design (*.json)")
        if path:
            with open(path, "w") as fh:
                json.dump(self.design.to_dict(), fh, indent=1)
