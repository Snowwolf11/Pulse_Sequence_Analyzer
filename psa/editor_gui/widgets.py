"""Shared widgets of the editor: themed pyqtgraph plots, draggable keyframes, matplotlib 3D view."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import Signal, Qt
from PySide6.QtWidgets import QVBoxLayout, QWidget

PANEL_BG = "#111827"
TEXT_FG = "#e5e7eb"
MUTED_FG = "#94a3b8"
ACCENT = "#38bdf8"
ACCENT_2 = "#f59e0b"
GOOD = "#22c55e"
BAD = "#ef4444"

pg.setConfigOptions(antialias=True, background=PANEL_BG, foreground=TEXT_FG)


def make_plot(y_label, x_label="time [µs]", y_units=None):
    w = pg.PlotWidget()
    w.showGrid(x=True, y=True, alpha=0.25)
    w.setLabel("left", y_label, units=y_units)
    w.setLabel("bottom", x_label)
    w.setMenuEnabled(False)
    return w


class KeyframeItem(pg.GraphItem):
    """Keyframe points of one channel: drag to move (snapped to sub-pulse starts), click to select,
    right-click to delete. Positions are given in plot coordinates (x in µs)."""

    sigDragged = Signal(int, float, float)     # position in list, x, y
    sigDragFinished = Signal()
    sigSelected = Signal(int)
    sigDeleteRequested = Signal(int)

    def __init__(self):
        super().__init__()
        self._drag = None
        self.selected = None

    def set_points(self, x, y, selected=None):
        self.selected = selected
        n = len(x)
        brushes = [pg.mkBrush(ACCENT_2 if i == selected else ACCENT) for i in range(n)]
        self.setData(pos=np.column_stack([x, y]) if n else np.zeros((0, 2)), size=11, symbol="o",
                     symbolBrush=brushes, symbolPen=pg.mkPen("w", width=1), pxMode=True, data=np.arange(n))

    def _point_at(self, pos):
        pts = self.scatter.pointsAt(pos)
        return int(pts[0].data()) if len(pts) else None

    def mouseDragEvent(self, ev):
        if ev.button() != Qt.MouseButton.LeftButton:
            ev.ignore()
            return
        if ev.isStart():
            i = self._point_at(ev.buttonDownPos())
            if i is None:
                ev.ignore()
                return
            self._drag = i
            self.sigSelected.emit(i)
        elif self._drag is None:
            ev.ignore()
            return
        p = ev.pos()
        self.sigDragged.emit(self._drag, float(p.x()), float(p.y()))
        if ev.isFinish():
            self._drag = None
            self.sigDragFinished.emit()
        ev.accept()

    def mouseClickEvent(self, ev):
        i = self._point_at(ev.pos())
        if i is None:
            ev.ignore()
            return
        if ev.button() == Qt.MouseButton.RightButton:
            self.sigDeleteRequested.emit(i)
        else:
            self.sigSelected.emit(i)
        ev.accept()


class Mpl3D(QWidget):
    """Small matplotlib 3D view."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(4, 3.4))
        self.figure.patch.set_facecolor(PANEL_BG)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111, projection="3d")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)
        self._style()

    def _style(self):
        ax = self.ax
        ax.set_facecolor(PANEL_BG)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_pane_color((0.07, 0.09, 0.15, 1.0))
            axis.label.set_color(TEXT_FG)
        ax.tick_params(colors=MUTED_FG, labelsize=7)

    def plot_curve(self, curve, title="", scale=1e6, unit="µs"):
        ax = self.ax
        ax.clear()
        self._style()
        if curve is not None and len(curve):
            c = np.asarray(curve) * scale
            ax.plot(c[:, 0], c[:, 1], c[:, 2], color=ACCENT, lw=1.2)
            ax.scatter(*c[0], color=GOOD, s=18)
            ax.scatter(*c[-1], color=BAD, s=18)
            span = np.ptp(c, axis=0).max() / 2 or 1.0
            mid = (c.max(axis=0) + c.min(axis=0)) / 2
            ax.set_xlim(mid[0] - span, mid[0] + span)
            ax.set_ylim(mid[1] - span, mid[1] + span)
            ax.set_zlim(mid[2] - span, mid[2] + span)
        ax.set_xlabel(f"x [{unit}]", fontsize=8)
        ax.set_ylabel(f"y [{unit}]", fontsize=8)
        ax.set_zlabel(f"z [{unit}]", fontsize=8)
        ax.set_title(title, color=TEXT_FG, fontsize=9)
        self.canvas.draw_idle()
