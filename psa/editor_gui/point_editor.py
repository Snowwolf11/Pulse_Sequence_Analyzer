"""Visual editor of the GRAPE point selection on the offset × B1 plane (separate, non-modal window).

Mouse: click on empty space adds a point (with the weight of the weight field; mirrored to −offset if
"± symmetric" is on), click on a point selects it (Shift / Ctrl: add to the selection), drag a point moves the
selection, drag on empty space selects in a rectangle (or fills it with a grid, see "drag on empty space"),
right click on a point deletes it, Shift + wheel over a point changes its weight (×1.25 per step), wheel zooms,
right drag scales the axes. Keys: Del / Backspace delete, Ctrl+A select all, Ctrl+Z / Ctrl+Shift+Z undo / redo.
Lock to grid snaps points to the steps of Settings → Solver (default 0.1 kHz and 0.5 %).
Background (optional): log10 infidelity of the tab's current pulse (start or result) for the gate goal.
Marker size grows with the weight; the table gives the exact values (editable).
"""
from __future__ import annotations

import math

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QItemSelectionModel, QRectF, Qt, Signal
from PySide6.QtGui import QKeySequence, QPen, QShortcut
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QGraphicsRectItem, QHBoxLayout, QHeaderView, QLabel, QPushButton,
                               QSpinBox, QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from psa.design.points import DEFAULT_GRID, load_points, save_points, snap
from psa.editor_gui.widgets import ACCENT, BAD, MUTED_FG, lock_heights

PICK_PX = 9                       # pick radius around a point [pixels]
WHEEL_FACTOR = 1.25
BG_LEVELS = (-8.0, 0.0)


def marker_size(w):
    return float(np.clip(10 + 4 * math.log2(max(w, 1e-6)), 6, 28))


class _PlaneView(pg.ViewBox):
    """ViewBox that hands left clicks / drags and Shift+wheel to the editor (pan / zoom otherwise)."""

    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        self.setMouseEnabled(True, True)

    def mouseClickEvent(self, ev):
        p = self.mapSceneToView(ev.scenePos())
        if ev.button() == Qt.MouseButton.LeftButton:
            ev.accept()
            self.editor.click(p.x(), p.y(), ev.modifiers())
        elif ev.button() == Qt.MouseButton.RightButton and self.editor.point_at(p.x(), p.y()) is not None:
            ev.accept()
            self.editor.delete_at(p.x(), p.y())
        else:
            super().mouseClickEvent(ev)

    def mouseDragEvent(self, ev, axis=None):
        if ev.button() != Qt.MouseButton.LeftButton:
            super().mouseDragEvent(ev, axis)
            return
        ev.accept()
        p0, p = self.mapSceneToView(ev.buttonDownScenePos()), self.mapSceneToView(ev.scenePos())
        self.editor.drag(p0.x(), p0.y(), p.x(), p.y(), ev.isStart(), ev.isFinish(), ev.modifiers())

    def wheelEvent(self, ev, axis=None):
        if ev.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            p = self.mapSceneToView(ev.scenePos())
            i = self.editor.point_at(p.x(), p.y())
            if i is not None:
                ev.accept()
                self.editor.scale_weight(i, WHEEL_FACTOR if ev.delta() > 0 else 1 / WHEEL_FACTOR)
                return
        super().wheelEvent(ev, axis)


class PointEditor(QDialog):
    """Edits a point selection; ``applied`` carries [[offset_hz, scale, weight], ...]."""
    applied = Signal(object)

    def __init__(self, points=None, grid_provider=None, background_provider=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("GRAPE points: offset × B1")
        self.resize(1100, 640)
        self.grid_provider = grid_provider or (lambda: dict(DEFAULT_GRID))
        self.background_provider = background_provider      # callable() -> (pulse, GrapeConfig) or None
        self.pts = []                                       # [offset_khz, b1_percent, weight]
        self.sel = set()
        self._undo, self._redo = [], []
        self._drag = None
        self._updating_table = False

        outer = QHBoxLayout(self)
        split = QSplitter(Qt.Orientation.Horizontal); outer.addWidget(split)
        # ---- plane
        self.view = _PlaneView(self)
        self.plot = pg.PlotWidget(viewBox=self.view)
        self.plot.setLabel("bottom", "offset [kHz]"); self.plot.setLabel("left", "B1 [%]")
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        self.bg_item = pg.ImageItem(); self.bg_item.setZValue(-10); self.bg_item.setVisible(False)
        self.plot.addItem(self.bg_item)
        try:
            cmap = pg.colormap.get("jet", source="matplotlib")          # the app's colormap
        except Exception:                                   # matplotlib missing: a pyqtgraph map
            cmap = pg.colormap.get("viridis")
        self.colorbar = pg.ColorBarItem(values=BG_LEVELS, colorMap=cmap, interactive=False, width=12)
        self.colorbar.setImageItem(self.bg_item, insert_in=self.plot.getPlotItem())
        self.colorbar.setVisible(False)
        self.scatter = pg.ScatterPlotItem(pxMode=True)
        self.plot.addItem(self.scatter)
        self.band = QGraphicsRectItem(); self.band.setPen(QPen(pg.mkColor(ACCENT), 0)); self.band.setVisible(False)
        self.band.setBrush(pg.mkBrush(56, 189, 248, 40))
        self.view.addItem(self.band, ignoreBounds=True)
        self.hover = QLabel(""); self.hover.setStyleSheet(f"color: {MUTED_FG};")
        self.plot.scene().sigMouseMoved.connect(self._mouse_moved)
        pw = QWidget(); pl = QVBoxLayout(pw); pl.setContentsMargins(0, 0, 0, 0)
        pl.addWidget(self.plot, 1); pl.addWidget(self.hover)
        split.addWidget(pw)

        # ---- controls
        side = QWidget(); sl = QVBoxLayout(side); sl.setContentsMargins(4, 0, 0, 0)
        f = QFormLayout()
        self.weight_spin = QDoubleSpinBox(); self.weight_spin.setRange(1e-3, 1e3); self.weight_spin.setDecimals(3)
        self.weight_spin.setValue(1.0); self.weight_spin.setSingleStep(0.5)
        self.weight_spin.setToolTip("weight of new points; 'Set' applies it to the selected points")
        set_w = QPushButton("Set for selected"); set_w.clicked.connect(lambda: self.set_weight_selected(self.weight_spin.value()))
        wrow = QHBoxLayout(); wrow.addWidget(self.weight_spin); wrow.addWidget(set_w)
        f.addRow("weight:", wrow)
        self.drag_box = QComboBox(); self.drag_box.addItems(["select (rectangle)", "fill rectangle with a grid"])
        f.addRow("drag on empty space:", self.drag_box)
        frow = QHBoxLayout()
        self.fill_off = QSpinBox(); self.fill_off.setRange(1, 1000); self.fill_off.setValue(11); self.fill_off.setSuffix(" offsets")
        self.fill_b1 = QSpinBox(); self.fill_b1.setRange(1, 1000); self.fill_b1.setValue(3); self.fill_b1.setSuffix(" B1")
        frow.addWidget(self.fill_off); frow.addWidget(QLabel("×")); frow.addWidget(self.fill_b1)
        f.addRow("fill grid:", frow)
        self.mirror_check = QCheckBox("± symmetric (new points also at −offset)")
        f.addRow(self.mirror_check)
        self.lock_check = QCheckBox("lock to grid")
        self.lock_check.toggled.connect(self._lock_toggled)
        self.lock_info = QLabel(""); self.lock_info.setStyleSheet(f"color: {MUTED_FG};")
        lrow = QHBoxLayout(); lrow.addWidget(self.lock_check); lrow.addWidget(self.lock_info); lrow.addStretch()
        f.addRow(lrow)
        self.bg_check = QCheckBox("background: log10(1 − F) of the current pulse")
        self.bg_check.setToolTip("log10(1 − F) of the tab's current pulse (result, else start) for the gate goal, over "
                                 "the visible area; 'Update' recomputes it for the current view (jet: −8 blue … 0 red)")
        self.bg_check.toggled.connect(self._bg_toggled)
        bg_upd = QPushButton("Update"); bg_upd.clicked.connect(self.update_background)
        brow = QHBoxLayout(); brow.addWidget(self.bg_check); brow.addWidget(bg_upd); brow.addStretch()
        f.addRow(brow)
        sl.addLayout(f)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["offset [kHz]", "B1 [%]", "weight"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.cellChanged.connect(self._cell_changed)
        self.table.itemSelectionChanged.connect(self._table_selection)
        sl.addWidget(self.table, 1)
        r1 = QHBoxLayout()
        for text, fn in (("Select all", self.select_all), ("Delete selected", self.delete_selected),
                         ("Undo", self.undo), ("Redo", self.redo)):
            b = QPushButton(text); b.clicked.connect(fn); r1.addWidget(b)
        sl.addLayout(r1)
        r2 = QHBoxLayout()
        for text, fn in (("Load…", self._load_dialog), ("Save…", self._save_dialog), ("Fit view", self.fit_view)):
            b = QPushButton(text); b.clicked.connect(fn); r2.addWidget(b)
        r2.addStretch()
        sl.addLayout(r2)
        self.info = QLabel(""); self.info.setWordWrap(True)
        sl.addWidget(self.info)
        r3 = QHBoxLayout(); r3.addStretch()
        self.apply_button = QPushButton("Apply"); self.apply_button.clicked.connect(self.apply)
        close = QPushButton("Close"); close.clicked.connect(self.hide)
        r3.addWidget(self.apply_button); r3.addWidget(close)
        sl.addLayout(r3)
        split.addWidget(side)
        split.setStretchFactor(0, 3); split.setStretchFactor(1, 2)

        for keys, fn in ((QKeySequence.StandardKey.Delete, self.delete_selected), ("Backspace", self.delete_selected),
                         (QKeySequence.StandardKey.SelectAll, self.select_all), (QKeySequence.StandardKey.Undo, self.undo),
                         (QKeySequence.StandardKey.Redo, self.redo), ("Ctrl+Shift+Z", self.redo)):
            sc = QShortcut(QKeySequence(keys), self); sc.activated.connect(fn)
        lock_heights(self)
        self.set_points(points or [])
        self.fit_view()

    # ------------------------------------------------------------ model ---
    def set_points(self, points):
        """Load a selection given as [[offset_hz, scale, weight], ...] (clears undo)."""
        self.pts = [[p[0] / 1e3, p[1] * 100.0, p[2] if len(p) > 2 else 1.0] for p in points]
        self.sel = set()
        self._undo, self._redo = [], []
        self._refresh()

    def points(self):
        """The selection as [[offset_hz, scale, weight], ...] (sorted by offset, then B1)."""
        return [[o * 1e3, b / 100.0, w] for o, b, w in sorted(self.pts, key=lambda p: (p[0], p[1]))]

    def grid(self):
        g = self.grid_provider() or {}
        return {"offset_khz": float(g.get("offset_khz", DEFAULT_GRID["offset_khz"])),
                "b1_percent": float(g.get("b1_percent", DEFAULT_GRID["b1_percent"]))}

    def _snap(self, o, b):
        return snap(o, b, self.grid()) if self.lock_check.isChecked() else (o, b)

    def _push(self):
        self._undo.append(([list(p) for p in self.pts], set(self.sel)))
        del self._undo[:-200]
        self._redo = []

    def undo(self):
        if self._undo:
            self._redo.append(([list(p) for p in self.pts], set(self.sel)))
            self.pts, self.sel = self._undo.pop()
            self._refresh()

    def redo(self):
        if self._redo:
            self._undo.append(([list(p) for p in self.pts], set(self.sel)))
            self.pts, self.sel = self._redo.pop()
            self._refresh()

    def _new(self, o, b, w):
        """Add a point unless one is already there (returns its index)."""
        o, b = self._snap(o, b)
        if b <= 0:
            return None
        for i, p in enumerate(self.pts):
            if abs(p[0] - o) < 1e-9 and abs(p[1] - b) < 1e-9:
                return i
        self.pts.append([o, b, w])
        return len(self.pts) - 1

    def add_point(self, o, b, w=None, mirror=None):
        """Add a point at (offset kHz, B1 %) (and at −offset if mirrored); selects it."""
        w = self.weight_spin.value() if w is None else w
        mirror = self.mirror_check.isChecked() if mirror is None else mirror
        self._push()
        idx = [self._new(o, b, w)] + ([self._new(-o, b, w)] if mirror and abs(o) > 1e-12 else [])
        self.sel = {i for i in idx if i is not None}
        self._refresh()

    def point_at(self, x, y):
        """Index of the point under the view position (within PICK_PX pixels), or None."""
        if not self.pts:
            return None
        px, py = self._pixel_size()
        d = [math.hypot((p[0] - x) / px, (p[1] - y) / py) for p in self.pts]
        i = int(np.argmin(d))
        return i if d[i] <= PICK_PX else None

    def _pixel_size(self):
        r = self.view.viewRect()
        w, h = max(self.view.width(), 1.0), max(self.view.height(), 1.0)
        return r.width() / w, r.height() / h

    def select(self, indices, extend=False):
        self.sel = (self.sel if extend else set()) | set(indices)
        self._refresh()

    def select_all(self):
        self.select(range(len(self.pts)))

    def select_in_rect(self, o0, b0, o1, b1, extend=False):
        lo_o, hi_o, lo_b, hi_b = min(o0, o1), max(o0, o1), min(b0, b1), max(b0, b1)
        self.select([i for i, p in enumerate(self.pts) if lo_o <= p[0] <= hi_o and lo_b <= p[1] <= hi_b], extend)

    def fill_rect(self, o0, b0, o1, b1, n_off=None, n_b1=None):
        """Fill the rectangle with an n_off × n_b1 grid of points (weight of the weight field)."""
        n_off = self.fill_off.value() if n_off is None else n_off
        n_b1 = self.fill_b1.value() if n_b1 is None else n_b1
        offs = np.linspace(min(o0, o1), max(o0, o1), n_off) if n_off > 1 else [0.5 * (o0 + o1)]
        b1s = np.linspace(min(b0, b1), max(b0, b1), n_b1) if n_b1 > 1 else [0.5 * (b0 + b1)]
        self._push()
        w, idx = self.weight_spin.value(), set()
        for o in offs:
            for b in b1s:
                idx.add(self._new(float(o), float(b), w))
                if self.mirror_check.isChecked() and abs(o) > 1e-12:
                    idx.add(self._new(-float(o), float(b), w))
        self.sel = {i for i in idx if i is not None}
        self._refresh()

    def move_selected(self, do, db, push=True):
        if not self.sel:
            return
        if push:
            self._push()
        for i in self.sel:
            o, b = self._snap(self.pts[i][0] + do, self.pts[i][1] + db)
            if b > 0:
                self.pts[i][0], self.pts[i][1] = o, b
        self._refresh()

    def delete_selected(self):
        if self.sel:
            self._push()
            self.pts = [p for i, p in enumerate(self.pts) if i not in self.sel]
            self.sel = set()
            self._refresh()

    def delete_at(self, x, y):
        i = self.point_at(x, y)
        if i is not None:
            self.sel = {i}
            self.delete_selected()

    def set_weight_selected(self, w):
        if self.sel and w > 0:
            self._push()
            for i in self.sel:
                self.pts[i][2] = float(w)
            self._refresh()

    def scale_weight(self, i, factor):
        self._push()
        self.pts[i][2] = float(np.clip(self.pts[i][2] * factor, 1e-3, 1e3))
        self.sel = {i}
        self._refresh()

    def _lock_toggled(self, on):
        g = self.grid()
        self.lock_info.setText(f"({g['offset_khz']:g} kHz, {g['b1_percent']:g} %; Settings → Solver)")
        if on and self.pts:                   # snap the existing points; merge those that land on the same spot
            self._push()
            merged = {}
            for i, (o, b, w) in enumerate(self.pts):
                key = snap(o, b, g)
                if key[1] > 0 and key not in merged:
                    merged[key] = [key[0], key[1], w]
            self.pts = list(merged.values())
            self.sel = set()
            self._refresh()

    # ------------------------------------------------------------ mouse ---
    def click(self, x, y, modifiers=Qt.KeyboardModifier.NoModifier):
        extend = bool(modifiers & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier))
        i = self.point_at(x, y)
        if i is None:
            self.add_point(x, y)
        elif extend and i in self.sel:
            self.select(self.sel - {i})
        else:
            self.select([i], extend)

    def drag(self, x0, y0, x, y, start, finish, modifiers=Qt.KeyboardModifier.NoModifier):
        if start:
            i = self.point_at(x0, y0)
            if i is not None:
                if i not in self.sel:
                    self.sel = {i}
                self._push()
                self._drag = {"kind": "move", "orig": {k: list(self.pts[k]) for k in self.sel}}
            else:
                self._drag = {"kind": "rect"}
        if self._drag is None:
            return
        if self._drag["kind"] == "move":
            for k, (o, b, _) in self._drag["orig"].items():
                o2, b2 = self._snap(o + x - x0, b + y - y0)
                if b2 > 0:
                    self.pts[k][0], self.pts[k][1] = o2, b2
            self._refresh(table=finish)
        else:
            self.band.setRect(QRectF(min(x0, x), min(y0, y), abs(x - x0), abs(y - y0)).normalized())
            self.band.setVisible(not finish)
            if finish:
                if self.drag_box.currentIndex() == 1:
                    self.fill_rect(x0, y0, x, y)
                else:
                    self.select_in_rect(x0, y0, x, y, bool(modifiers & (Qt.KeyboardModifier.ShiftModifier |
                                                                        Qt.KeyboardModifier.ControlModifier)))
        if finish:
            self._drag = None

    def _mouse_moved(self, pos):
        if self.plot.sceneBoundingRect().contains(pos):
            p = self.view.mapSceneToView(pos)
            i = self.point_at(p.x(), p.y())
            if i is not None:
                o, b, w = self.pts[i]
                self.hover.setText(f"point: {o:g} kHz, {b:g} %, weight {w:g}")
            else:
                o, b = self._snap(p.x(), p.y())
                self.hover.setText(f"{o:.3f} kHz, {b:.2f} %")

    # ------------------------------------------------------------- view ---
    def _refresh(self, table=True):
        spots = []
        for i, (o, b, w) in enumerate(self.pts):
            s = i in self.sel
            spots.append({"pos": (o, b), "size": marker_size(w),
                          "brush": pg.mkBrush(ACCENT if s else "#f5f5f4"),
                          "pen": pg.mkPen("#f59e0b" if s else "#1f2937", width=2 if s else 1)})
        self.scatter.setData(spots)
        if table:
            self._fill_table()
        n, ws = len(self.pts), [p[2] for p in self.pts]
        wtxt = "" if not n else (", weights 1" if min(ws) == max(ws) == 1 else f", weights {min(ws):g} … {max(ws):g}")
        self.info.setText(f"{n} point{'s' * (n != 1)}{wtxt}, {len(self.sel)} selected")
        self.info.setStyleSheet("" if n else f"color: {BAD};")
        self.apply_button.setEnabled(n > 0)

    def _fill_table(self):
        self._updating_table = True
        try:
            self.table.setRowCount(len(self.pts))
            for r, p in enumerate(self.pts):
                for c, v in enumerate(p):
                    self.table.setItem(r, c, QTableWidgetItem(f"{v:g}"))
            self.table.clearSelection()
            sm, model = self.table.selectionModel(), self.table.model()
            flags = QItemSelectionModel.SelectionFlag.Select | QItemSelectionModel.SelectionFlag.Rows
            for r in sorted(self.sel):
                sm.select(model.index(r, 0), flags)
        finally:
            self._updating_table = False

    def _cell_changed(self, r, c):
        if self._updating_table or r >= len(self.pts):
            return
        try:
            v = float(self.table.item(r, c).text().replace(",", "."))
            if c > 0 and v <= 0:
                raise ValueError
        except ValueError:
            self._fill_table()
            return
        self._push()
        self.pts[r][c] = v
        self._refresh()

    def _table_selection(self):
        if self._updating_table:
            return
        self.sel = {i.row() for i in self.table.selectionModel().selectedRows()}
        self._refresh(table=False)

    def fit_view(self):
        if self.pts:
            o = [p[0] for p in self.pts]; b = [p[1] for p in self.pts]
            po, pb = max(1.0, 0.1 * (max(o) - min(o))), max(2.0, 0.1 * (max(b) - min(b)))
            self.view.setRange(xRange=(min(o) - po, max(o) + po), yRange=(min(b) - pb, max(b) + pb), padding=0)
        else:
            self.view.setRange(xRange=(-12, 12), yRange=(80, 120), padding=0)

    def _bg_toggled(self, on):
        if on:
            self.update_background()
        self.bg_item.setVisible(on and self.bg_item.image is not None)
        self.colorbar.setVisible(self.bg_item.isVisible())

    def background_map(self, nx=161, ny=61):
        """(offsets_khz, b1_percent, log10 infidelity [ny, nx]) of the current pulse over the visible area,
        or None if there is no pulse / gate goal."""
        from psa.design.grape import band_infidelity
        src = self.background_provider() if self.background_provider else None
        if src is None:
            return None
        ps, cfg = src
        r = self.view.viewRect()
        offs = np.linspace(r.left(), r.right(), nx)
        b1 = np.linspace(max(r.top(), 1e-3), max(r.bottom(), 2e-3), ny)
        z = np.array([band_infidelity(ps, cfg, offs * 1e3, scale=b / 100.0) for b in b1])
        return offs, b1, np.log10(np.maximum(z, 1e-16))

    def update_background(self):
        try:
            m = self.background_map()
        except Exception as e:                      # e.g. invalid gate inputs
            self.hover.setText(f"background: {e}")
            m = None
        if m is None:
            self.bg_item.setVisible(False)
            self.colorbar.setVisible(False)
            if self.bg_check.isChecked():
                self.hover.setText("background: needs a pulse and a gate goal (UR or PP)")
            return
        offs, b1, z = m
        self.bg_item.setImage(z.T, levels=BG_LEVELS, autoLevels=False)
        do, db = offs[1] - offs[0], b1[1] - b1[0]
        self.bg_item.setRect(QRectF(offs[0] - do / 2, b1[0] - db / 2, do * len(offs), db * len(b1)))
        self.bg_item.setVisible(self.bg_check.isChecked())
        self.colorbar.setVisible(self.bg_item.isVisible())

    # ------------------------------------------------------------ files ---
    def save(self, path):
        save_points(path, self.points(), self.grid())

    def load(self, path):
        pts, _ = load_points(path)
        self._push()
        self.pts = [[p[0] / 1e3, p[1] * 100.0, p[2]] for p in pts]
        self.sel = set()
        self._refresh()
        self.fit_view()

    def _save_dialog(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save point selection", "points.json",
                                              "Point selection (*.json);;CSV (*.csv)")
        if path:
            try:
                self.save(path)
            except OSError as e:
                self.info.setText(f"could not save: {e}"); self.info.setStyleSheet(f"color: {BAD};")

    def _load_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load point selection", "", "Point selection (*.json *.csv)")
        if path:
            try:
                self.load(path)
            except (OSError, ValueError, KeyError) as e:
                self.info.setText(f"could not load: {e}"); self.info.setStyleSheet(f"color: {BAD};")

    def apply(self):
        if self.pts:
            self.applied.emit(self.points())
