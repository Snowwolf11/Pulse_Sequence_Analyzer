"""Animation tab of the analyzer: the spin dynamics of the main window's pulse in the rotating frame.

Four synchronized 3D panels (each can be hidden; the Bloch-sphere views share one camera):
    Ensemble         every spin of the ensemble (offset x B1) on the Bloch sphere: colour = offset (blue ->
                     red), size = B1, fading trails; the selected spin is ringed; optional target
    Effective field  the selected spin: magnetization (cyan), effective field of the running sub-pulse
                     (orange) and the precession circle around it
    Rotation         the selected spin's rotation R(t): the rotated frame x / y / z (with the final frame as
                     ghosts), the DROPS droplets of the propagator U (rank 0 and rank 1, colour = phase), or the
                     Q-bead of the magnetization (sphere coloured red / green by m . r)
    Error curve      the selected spin's error curve growing with time (green -> red), current tangent
Below: the pulse (amplitude, phase) with a draggable time cursor, playback controls and the settings.
The dynamics (psa.design.animation) are computed once per pulse / setting; playback only draws. Export GIF
writes the visible panels side by side. Needs PyOpenGL (otherwise a note is shown).
"""
from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QRectF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QImage, QLinearGradient, QPainter
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFrame, QGridLayout, QHBoxLayout,
                               QLabel, QLineEdit, QPushButton, QSlider, QSpinBox, QVBoxLayout, QWidget)

from psa.design.animation import (Dynamics, droplet, magnetization_droplet, phase_colors, propagator_droplets,
                                  rotate, unit_sphere)
from psa.design.points import parse_values

try:
    import pyqtgraph.opengl as gl
    from OpenGL import GL  # noqa: F401
    HAVE_GL = True
except Exception:                                    # PyOpenGL missing
    gl, HAVE_GL = None, False

BG = (14, 18, 30, 255)
SPHERE = (0.55, 0.65, 0.90, 0.07)
WIRE = (0.65, 0.72, 0.85, 0.35)
CYAN, ORANGE, WHITE = (0.30, 0.85, 1.0, 1.0), (1.0, 0.62, 0.15, 1.0), (1.0, 1.0, 1.0, 1.0)
FRAME_COLORS = ((0.95, 0.30, 0.30, 1.0), (0.35, 0.85, 0.40, 1.0), (0.35, 0.55, 1.0, 1.0))


def opengl_works():
    """True if an OpenGL context can be created here (not e.g. on Qt's offscreen platform or without drivers):
    showing a GL window without one crashes inside Qt."""
    if not HAVE_GL:
        return False
    try:
        from PySide6.QtGui import QGuiApplication, QOpenGLContext
        if QGuiApplication.platformName() in ("offscreen", "minimal"):
            return False
        ctx = QOpenGLContext()
        return bool(ctx.create())
    except Exception:
        return False


def offset_colors(values, lo, hi, alpha=1.0):
    """Blue (lo) -> white-ish -> red (hi), RGBA (n, 4)."""
    import matplotlib
    cmap = matplotlib.colormaps["coolwarm"] if hasattr(matplotlib, "colormaps") else matplotlib.cm.get_cmap("coolwarm")
    x = (np.asarray(values, float) - lo) / (hi - lo) if hi > lo else np.full(len(values), 0.5)
    c = cmap(np.clip(x, 0, 1))
    c[:, 3] = alpha
    return c


def time_colors(n, alpha=1.0):
    """Green (start) -> yellow -> red (end), RGBA (n, 4)."""
    x = np.linspace(0, 1, max(n, 1))[:, None]
    a, b, c = np.array([0.25, 0.85, 0.35]), np.array([0.95, 0.85, 0.25]), np.array([0.95, 0.30, 0.25])
    rgb = np.where(x < 0.5, a + (b - a) * (2 * x), b + (c - b) * (2 * x - 1))
    return np.column_stack([rgb, np.full(len(rgb), alpha)])


class ColorBar(QWidget):
    """A small horizontal colour bar with labels at both ends."""

    def __init__(self, lo_text="", hi_text="", parent=None):
        super().__init__(parent)
        self.lo_text, self.hi_text = lo_text, hi_text
        self.setFixedHeight(16)

    def set_labels(self, lo, hi):
        self.lo_text, self.hi_text = lo, hi
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        w = self.width()
        bar = QRectF(60, 3, max(w - 120, 10), 10)
        g = QLinearGradient(bar.left(), 0, bar.right(), 0)
        for x, c in zip(np.linspace(0, 1, 9), offset_colors(np.linspace(0, 1, 9), 0, 1)):
            g.setColorAt(float(x), QColor.fromRgbF(*c[:3]))
        p.fillRect(bar, g)
        p.setPen(QColor(200, 205, 215))
        p.drawText(QRectF(0, 0, 58, 16), Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, self.lo_text)
        p.drawText(QRectF(bar.right() + 2, 0, 58, 16), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                   self.hi_text)


if HAVE_GL:
    from OpenGL.GL import glViewport

    from PySide6.QtGui import QSurfaceFormat
    from PySide6.QtOpenGL import QOpenGLWindow

    _G = gl.GLViewWidget
    import inspect as _inspect
    NEW_PG = "viewport" in _inspect.signature(_G.setProjection).parameters     # pyqtgraph >= 0.14

    class _MultiViewBase:
        """Several panels side by side in ONE OpenGL surface, each in its own viewport with its own camera and
        items (one composite per frame instead of one per panel). The cameras of the ``synced`` panels move
        together. Used with a QOpenGLWidget (default, reliable) or a native QOpenGLWindow (faster in big
        windows on macOS, experimental: PSA_GL_NATIVE=1)."""

        def _setup(self, keys, synced):
            self.keys = list(keys)
            self.panel_items = {k: [] for k in self.keys}
            self.cams = {k: dict(distance=2.9, elevation=22.0, azimuth=35.0) for k in self.keys}
            self.shown = {k: True for k in self.keys}
            self.synced = set(synced)
            self._active = None
            self.opts["fov"] = 50
            self.setBackgroundColor(BG)

        @staticmethod
        def _format():
            fmt = QSurfaceFormat(); fmt.setSamples(4); fmt.setDepthBufferSize(24)     # smooth lines and edges
            return fmt

        def add(self, key, item):
            self.addItem(item)                        # sets the item's view, initializes it if possible
            self.items.remove(item)
            self.panel_items[key].append(item)

        def gl_ok(self):
            """True if a working OpenGL context exists (without one every GL call would crash)."""
            ctx = self.context()
            return ctx is not None and ctx.isValid()

        def initializeGL(self):
            if not self.gl_ok():
                return
            if isinstance(self, gl.GLViewWidget):         # pyqtgraph's own setup (0.14: the default VAO)
                gl.GLViewWidget.initializeGL(self)
            for items in self.panel_items.values():
                for it in items:
                    if not it.isInitialized():
                        it.initialize()

        def visible_keys(self):
            return [k for k in self.keys if self.shown[k]]

        def panel_rect(self, i, n):
            w, h = self.deviceWidth(), self.deviceHeight()
            x0, x1 = int(round(i * w / n)), int(round((i + 1) * w / n))
            return x0, 0, x1 - x0, h

        def _use(self, key):
            c = self.cams[key]
            self.opts["distance"], self.opts["elevation"], self.opts["azimuth"] = (c["distance"], c["elevation"],
                                                                                   c["azimuth"])

        def paintGL(self):
            from OpenGL.GL import GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, glClear, glClearColor
            if not self.gl_ok():
                return
            keys = self.visible_keys()
            glViewport(0, 0, self.deviceWidth(), self.deviceHeight())
            glClearColor(*self.opts["bgcolor"]); glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            saved = self.items
            for i, k in enumerate(keys):
                vp = self.panel_rect(i, len(keys))
                self.opts["viewport"] = vp
                self._use(k)
                self.items = self.panel_items[k]
                glViewport(*vp)
                if NEW_PG:                                # pyqtgraph >= 0.14: logical pixels, matrix stacks
                    r = self.devicePixelRatioF() if hasattr(self, "devicePixelRatioF") else self.devicePixelRatio()
                    lvp = tuple(v / r for v in vp)
                    self.setProjection(lvp, lvp)
                else:                                     # older: the frustum of this panel (device pixels)
                    self.setProjection(region=vp)
                self.setModelview()
                glClear(GL_DEPTH_BUFFER_BIT)
                self.drawItemTree()
            self.items, self.opts["viewport"] = saved, None

        def _key_at(self, x):
            keys = self.visible_keys()
            if not keys:
                return None
            return keys[int(np.clip(x / max(self.width(), 1) * len(keys), 0, len(keys) - 1))]

        def _store(self, key):
            moved = {"distance": self.opts["distance"], "elevation": self.opts["elevation"],
                     "azimuth": self.opts["azimuth"]}
            for k in (self.synced if key in self.synced else {key}):
                self.cams[k].update(moved)

        def mousePressEvent(self, ev):
            self._active = self._key_at(ev.position().x())
            self.mousePos = ev.position()

        def mouseMoveEvent(self, ev):
            if self._active is None:
                return
            self._use(self._active)
            _G.mouseMoveEvent(self, ev)
            self._store(self._active)
            self.update()

        def wheelEvent(self, ev):
            k = self._key_at(ev.position().x())
            if k is None:
                return
            self._use(k)
            _G.wheelEvent(self, ev)
            self._store(k)
            self.update()

    class MultiViewWidget(_MultiViewBase, gl.GLViewWidget):
        def __init__(self, keys, synced=()):
            gl.GLViewWidget.__init__(self)
            self.setFormat(self._format())
            self._setup(keys, synced)
            self.container = self
            self.setMinimumSize(300, 220)

    class MultiViewWindow(_MultiViewBase, QOpenGLWindow):
        # camera and drawing code of pyqtgraph's GLViewWidget (works on any object with width / update /
        # isValid / device size)
        reset, addItem, removeItem, setBackgroundColor = _G.reset, _G.addItem, _G.removeItem, _G.setBackgroundColor
        getViewport, setProjection, projectionMatrix = _G.getViewport, _G.setProjection, _G.projectionMatrix
        setModelview, viewMatrix, drawItemTree = _G.setModelview, _G.viewMatrix, _G.drawItemTree
        setCameraPosition, cameraPosition, orbit, pan, pixelSize = (_G.setCameraPosition, _G.cameraPosition,
                                                                    _G.orbit, _G.pan, _G.pixelSize)
        setCameraParams, cameraParams = _G.setCameraParams, _G.cameraParams

        def __init__(self, keys, synced=()):
            QOpenGLWindow.__init__(self)
            self.setFormat(self._format())
            from pyqtgraph import Vector
            from PySide6.QtGui import QQuaternion
            self.opts = {"center": Vector(0, 0, 0), "rotation": QQuaternion(1, 0, 0, 0), "distance": 10.0, "fov": 50,
                         "elevation": 30, "azimuth": 45, "viewport": None, "rotationMethod": "euler"}
            self.reset()
            self.items = []
            self._setup(keys, synced)
            self.container = QWidget.createWindowContainer(self)
            self.container.setMinimumSize(300, 220)

        def deviceWidth(self):
            return int(self.width() * self.devicePixelRatio())

        def deviceHeight(self):
            return int(self.height() * self.devicePixelRatio())

        def isVisible(self):
            return self.container.isVisible()

    # the native window borrows pyqtgraph < 0.14 internals; experimental, only there
    MultiView = MultiViewWindow if os.environ.get("PSA_GL_NATIVE") == "1" and not NEW_PG else MultiViewWidget

    class Arrow:
        """Shaft (line) + head (cone) from a start point along a vector."""

        def __init__(self, view, color, width=3.0, head=0.12):
            self.head = head
            self.line = gl.GLLinePlotItem(pos=np.zeros((2, 3)), color=color, width=width, antialias=True)
            md = gl.MeshData.cylinder(rows=1, cols=20, radius=[head * 0.45, 0.0], length=head)
            self.cone = gl.GLMeshItem(meshdata=md, color=color, smooth=True, shader="shaded")
            view.addItem(self.line); view.addItem(self.cone)

        def set(self, start, vec):
            start, vec = np.asarray(start, float), np.asarray(vec, float)
            n = np.linalg.norm(vec)
            if n < 1e-9:
                self.line.setData(pos=np.array([start, start])); self.cone.setVisible(False); return
            d = vec / n
            base = start + vec - d * self.head
            self.line.setData(pos=np.array([start, base]))
            self.cone.setVisible(True)
            self.cone.resetTransform()
            z = np.array([0.0, 0.0, 1.0])
            ax = np.cross(z, d)
            ang = np.degrees(np.arctan2(np.linalg.norm(ax), d @ z))
            if np.linalg.norm(ax) < 1e-9:
                ax = np.array([1.0, 0.0, 0.0])
            self.cone.rotate(ang, *ax)
            self.cone.translate(*base)

        def setVisible(self, on):
            self.line.setVisible(on); self.cone.setVisible(on)


def mesh(verts, faces, colors=None):
    """MeshData with vertex normals computed vectorized (pyqtgraph loops over the vertices in Python: too slow to
    rebuild a droplet every frame)."""
    verts = np.ascontiguousarray(verts, dtype=np.float32)
    md = gl.MeshData(vertexes=verts, faces=faces, vertexColors=colors)
    v = verts[faces]
    fn = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])                    # area-weighted face normals
    vn = np.zeros_like(verts)
    for k in range(3):
        np.add.at(vn, faces[:, k], fn)
    n = np.linalg.norm(vn, axis=1, keepdims=True)
    vn = np.where(n > 1e-12, vn / np.where(n > 1e-12, n, 1.0), verts / np.maximum(np.linalg.norm(verts, axis=1,
                                                                                              keepdims=True), 1e-12))
    md._vertexNormals = vn.astype(np.float32)
    return md


def _circle(normal, center=(0, 0, 0), radius=1.0, n=96):
    nrm = np.asarray(normal, float); nrm = nrm / np.linalg.norm(nrm)
    a = np.cross(nrm, [1.0, 0, 0]) if abs(nrm[0]) < 0.9 else np.cross(nrm, [0, 1.0, 0])
    a /= np.linalg.norm(a); b = np.cross(nrm, a)
    t = np.linspace(0, 2 * np.pi, n)
    return np.asarray(center, float) + radius * (np.outer(np.cos(t), a) + np.outer(np.sin(t), b))


class Panel:
    """One panel of the MultiView: its key, title / footer labels (above / below the 3D view) and items."""

    def __init__(self, key, title, view):
        self.key, self.mv = key, view
        self.title = QLabel(title); self.title.setStyleSheet("font-weight: bold;")
        self.footer = QLabel(""); self.footer.setStyleSheet("color: #94a3b8;"); self.footer.setWordWrap(True)
        self.footer.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.box = QWidget(); bl = QVBoxLayout(self.box); bl.setContentsMargins(4, 0, 4, 0)
        bl.addWidget(self.footer)

    def addItem(self, item):
        self.mv.add(self.key, item)

    @property
    def view(self):
        return self                               # items are added to the panel (Arrow(panel.view, ...))

    def isVisible(self):
        return self.mv.shown[self.key] and self.mv.isVisible()

    def setVisible(self, on):
        self.mv.shown[self.key] = bool(on)
        self.title.setVisible(on); self.box.setVisible(on)
        self.mv.update()

    def add_bloch_sphere(self):
        dirs, faces = unit_sphere(24, 48)
        md = gl.MeshData(vertexes=dirs, faces=faces)
        self.addItem(gl.GLMeshItem(meshdata=md, color=SPHERE, smooth=True, shader="shaded", glOptions="translucent"))
        for nrm in ((0, 0, 1), (1, 0, 0), (0, 1, 0)):
            self.addItem(gl.GLLinePlotItem(pos=_circle(nrm), color=WIRE, width=1.0, antialias=True))
        for k, vec in enumerate(((1.3, 0, 0), (0, 1.3, 0), (0, 0, 1.3))):        # axes x / y / z: red / green / blue
            c = FRAME_COLORS[k]
            self.addItem(gl.GLLinePlotItem(pos=np.array([(0, 0, 0), vec]), color=(c[0], c[1], c[2], 0.55), width=1.5))


class AnimationTab(QWidget):
    """The animation tab (see the module doc). ``pulse_provider()`` -> PulseSequence or None."""

    def __init__(self, pulse_provider, parent=None):
        super().__init__(parent)
        self.pulse_provider = pulse_provider
        self.dyn = None
        self.tau = 0.0
        self._dirty = True
        outer = QVBoxLayout(self)
        if not HAVE_GL:
            msg = QLabel("The animation needs OpenGL for Python. Install it with\n\n    python -m pip install PyOpenGL\n\n"
                         "and restart the analyzer.")
            msg.setAlignment(Qt.AlignmentFlag.AlignCenter)
            outer.addWidget(msg)
            self.panels = {}
            return
        # ---- panels: one GL widget, titles above and footers below in matching columns
        keys = (("ensemble", "Ensemble · rotating frame"), ("field", "Effective field · rotating frame"),
                ("rotation", "Rotation · rotating frame"), ("curve", "Error curve · toggling frame"))
        self.mv = MultiView([k for k, _ in keys], synced=("ensemble", "field", "rotation"))
        self.mv.cams["curve"]["distance"] = 2.6
        self.panels = {k: Panel(k, t, self.mv) for k, t in keys}
        trow, frow = QHBoxLayout(), QHBoxLayout()
        for p in self.panels.values():
            trow.addWidget(p.title, 1); frow.addWidget(p.box, 1)
        self.panel_row = QWidget(); pl = QVBoxLayout(self.panel_row); pl.setContentsMargins(0, 0, 0, 0)
        pl.setSpacing(2)
        self.gl_ready = opengl_works()
        self.no_gl = QLabel("OpenGL is not available here (no driver / remote session): the 3D panels cannot be "
                            "drawn. The playback controls and the GIF export need it too.")
        self.no_gl.setWordWrap(True); self.no_gl.setStyleSheet("color: #f87171;")
        self.no_gl.setVisible(not self.gl_ready)
        pl.addLayout(trow); pl.addWidget(self.no_gl)
        if self.gl_ready:
            pl.addWidget(self.mv.container, 1)
        else:                                     # never shown: a GL window without a context crashes Qt
            self.mv.container.hide()
        pl.addLayout(frow)
        outer.addWidget(self.panel_row, 5)
        self._build_ensemble(); self._build_field(); self._build_rotation(); self._build_curve()
        # ---- pulse strip with a time cursor
        self.pulse_plot = pg.PlotWidget(); self.pulse_plot.setMaximumHeight(130)
        self.pulse_plot.setLabel("left", "amplitude [%]"); self.pulse_plot.setLabel("bottom", "time [µs]")
        self.pulse_plot.showGrid(x=True, y=True, alpha=0.2)
        self.amp_curve = self.pulse_plot.plot(pen=pg.mkPen((96, 200, 255), width=1.5), stepMode="center")
        self.cursor = pg.InfiniteLine(angle=90, movable=True, pen=pg.mkPen((255, 160, 40), width=2))
        self.cursor.sigPositionChanged.connect(self._cursor_moved)
        self.pulse_plot.addItem(self.cursor)
        outer.addWidget(self.pulse_plot, 1)
        # ---- controls
        ctl = QHBoxLayout()
        self.play_button = QPushButton("▶ Play"); self.play_button.clicked.connect(self.toggle_play)
        self.slider = QSlider(Qt.Orientation.Horizontal); self.slider.setRange(0, 1000)
        self.slider.valueChanged.connect(self._slider_moved)
        self.time_label = QLabel("t = 0 µs"); self.time_label.setMinimumWidth(170)
        self.speed_spin = QDoubleSpinBox(); self.speed_spin.setRange(0.5, 600); self.speed_spin.setValue(8.0)
        self.speed_spin.setSuffix(" s"); self.speed_spin.setPrefix("pulse plays in ")
        self.loop_check = QCheckBox("loop"); self.loop_check.setChecked(True)
        self.trail_spin = QSpinBox(); self.trail_spin.setRange(0, 100); self.trail_spin.setValue(25)
        self.trail_spin.setPrefix("trail "); self.trail_spin.setSuffix(" % of the pulse")
        self.trail_spin.valueChanged.connect(lambda *a: self.draw())
        for w in (self.play_button, self.slider, self.time_label, self.speed_spin, self.loop_check, self.trail_spin):
            ctl.addWidget(w, 1 if w is self.slider else 0)
        outer.addLayout(ctl)
        # ---- settings
        g = QGridLayout()
        self.off_edit = QLineEdit("-10:10:11"); self.off_edit.setToolTip("offsets of the ensemble in kHz (a:b:n ranges)")
        self.b1_edit = QLineEdit("100"); self.b1_edit.setToolTip("B1 values of the ensemble in % (every offset with every value)")
        self.m0_spins = [self._spin(v, -1, 1, 2) for v in (0.0, 0.0, 1.0)]
        self.sel_off = self._spin(0.0, -1e4, 1e4, 2, " kHz"); self.sel_b1 = self._spin(100.0, 1, 1000, 1, " %")
        self.target_check = QCheckBox("target:")
        self.target_spins = [self._spin(v, -1, 1, 2) for v in (1.0, 0.0, 0.0)]
        self.rot_mode = QComboBox(); self.rot_mode.addItems(["frame x, y, z", "droplets of the rotation U",
                                                             "Q-bead of the magnetization"])
        self.rot_mode.currentIndexChanged.connect(self._rot_mode_changed)
        self.show_checks = {}
        g.addWidget(QLabel("ensemble offsets [kHz]:"), 0, 0); g.addWidget(self.off_edit, 0, 1)
        g.addWidget(QLabel("B1 [%]:"), 0, 2); g.addWidget(self.b1_edit, 0, 3)
        m0_label = QLabel("start m0:")
        m0_tip = ("magnetization of every spin at t = 0 (normalized), e.g. z for excitation, x to see what the pulse "
                  "does to transverse magnetization")
        m0_label.setToolTip(m0_tip)
        for sp in self.m0_spins:
            sp.setToolTip(m0_tip)
        g.addWidget(m0_label, 0, 4)
        m0row = QHBoxLayout(); [m0row.addWidget(s) for s in self.m0_spins]; g.addLayout(m0row, 0, 5)
        g.addWidget(QLabel("selected spin:"), 1, 0)
        selrow = QHBoxLayout(); selrow.addWidget(self.sel_off); selrow.addWidget(self.sel_b1); g.addLayout(selrow, 1, 1)
        g.addWidget(self.target_check, 1, 2)
        trow = QHBoxLayout(); [trow.addWidget(s) for s in self.target_spins]; g.addLayout(trow, 1, 3)
        g.addWidget(QLabel("rotation panel:"), 1, 4); g.addWidget(self.rot_mode, 1, 5)
        srow = QHBoxLayout(); srow.addWidget(QLabel("show:"))
        for k, p in self.panels.items():
            c = QCheckBox(p.title.text().split(" · ")[0]); c.setChecked(True)
            c.toggled.connect(lambda on, p=p: (p.setVisible(on), self._fit_panels(), self.draw()))
            self.show_checks[k] = c; srow.addWidget(c)
        srow.addStretch()
        self.export_fps = QSpinBox(); self.export_fps.setRange(5, 60); self.export_fps.setValue(20)
        self.export_fps.setSuffix(" fps")
        self.export_button = QPushButton("Export GIF…"); self.export_button.clicked.connect(self.export_gif_dialog)
        self.export_button.setToolTip("the visible panels side by side, one loop of the pulse in the playback time")
        srow.addWidget(self.export_fps); srow.addWidget(self.export_button)
        g.addLayout(srow, 2, 0, 1, 6)
        self.error = QLabel(""); self.error.setStyleSheet("color: #f87171;")
        g.addWidget(self.error, 3, 0, 1, 6)
        outer.addLayout(g)
        for w in [self.off_edit, self.b1_edit]:
            w.editingFinished.connect(self.recompute)
        for s in self.m0_spins + [self.sel_off, self.sel_b1]:
            s.valueChanged.connect(self._schedule)
        for s in self.target_spins:
            s.valueChanged.connect(lambda *a: self.draw())
        self.target_check.toggled.connect(lambda *a: self.draw())
        self._recompute_timer = QTimer(self); self._recompute_timer.setSingleShot(True); self._recompute_timer.setInterval(250)
        self._recompute_timer.timeout.connect(self.recompute)
        self.timer = QTimer(self); self.timer.setInterval(16); self.timer.timeout.connect(self._tick)
        self.timer.setTimerType(Qt.TimerType.PreciseTimer)
        self._last = None
        self._rot_mode_changed()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._fit_panels()

    def _fit_panels(self):
        """Keep the panels about square: the row is at most 1.25 x the width of one panel high."""
        if not self.panels:
            return
        n = max(1, len(self.mv.visible_keys()))
        self.panel_row.setMaximumHeight(int(1.25 * self.width() / n) + 60)

    @staticmethod
    def _spin(v, lo, hi, dec, suffix=""):
        s = QDoubleSpinBox(); s.setRange(lo, hi); s.setDecimals(dec); s.setValue(v); s.setSingleStep(0.1 if hi <= 1 else 1)
        if suffix:
            s.setSuffix(suffix)
        return s

    # ------------------------------------------------------------ build ---
    def _build_ensemble(self):
        p = self.panels["ensemble"]; p.add_bloch_sphere(); v = p.view
        self.ens_trails = gl.GLLinePlotItem(pos=np.zeros((2, 3)), mode="lines", width=2.0, antialias=True)
        self.ens_dots = gl.GLScatterPlotItem(pos=np.zeros((1, 3)), size=10, pxMode=True)
        self.ens_dots.setGLOptions("translucent")
        self.ens_sel = gl.GLScatterPlotItem(pos=np.zeros((1, 3)), size=30, color=(1, 1, 1, 0.30), pxMode=True)
        self.ens_target = gl.GLLinePlotItem(pos=np.zeros((2, 3)), color=(1, 0.9, 0.3, 1), width=2.5, antialias=True)
        for it in (self.ens_trails, self.ens_sel, self.ens_dots, self.ens_target):
            v.addItem(it)
        self.colorbar = ColorBar()
        p.box.layout().insertWidget(0, self.colorbar)

    def _build_field(self):
        p = self.panels["field"]; p.add_bloch_sphere(); v = p.view
        self.fld_circle = gl.GLLinePlotItem(pos=np.zeros((2, 3)), color=(1.0, 0.62, 0.15, 0.55), width=1.5, antialias=True)
        self.fld_trail = gl.GLLinePlotItem(pos=np.zeros((2, 3)), width=2.5, antialias=True)
        v.addItem(self.fld_circle); v.addItem(self.fld_trail)
        self.fld_b = Arrow(v, ORANGE, 3.5)
        self.fld_m = Arrow(v, CYAN, 3.5)

    def _build_rotation(self):
        p = self.panels["rotation"]; v = p.view
        self.rot_frame_items = []
        self.rot_arrows = [Arrow(v, c, 3.5) for c in FRAME_COLORS]
        self.rot_ghosts = [gl.GLLinePlotItem(pos=np.zeros((2, 3)), color=(c[0], c[1], c[2], 0.30), width=2.0)
                           for c in FRAME_COLORS]
        self.rot_tips = [gl.GLLinePlotItem(pos=np.zeros((2, 3)), color=(c[0], c[1], c[2], 0.45), width=1.5,
                                           antialias=True) for c in FRAME_COLORS]
        sphere = gl.GLMeshItem(meshdata=gl.MeshData(*unit_sphere(24, 48)), color=SPHERE, smooth=True,
                               shader="shaded", glOptions="translucent")
        for it in [sphere] + self.rot_ghosts + self.rot_tips:
            v.addItem(it)
        for k, vec in enumerate(((1.3, 0, 0), (0, 1.3, 0), (0, 0, 1.3))):
            c = FRAME_COLORS[k]
            ax = gl.GLLinePlotItem(pos=np.array([(0, 0, 0), vec]), color=(c[0], c[1], c[2], 0.35), width=1.0)
            v.addItem(ax); self.rot_frame_items.append(ax)
        self.rot_frame_items = [sphere] + self.rot_ghosts + self.rot_tips + self.rot_frame_items
        self.drop_dirs, self.drop_faces = unit_sphere(28, 56)
        self.drop0 = gl.GLMeshItem(meshdata=gl.MeshData(vertexes=self.drop_dirs, faces=self.drop_faces), smooth=True,
                                   shader="shaded")
        self.drop1 = gl.GLMeshItem(meshdata=gl.MeshData(vertexes=self.drop_dirs, faces=self.drop_faces), smooth=True,
                                   shader="shaded")
        self.drop_labels = []
        self.bead = gl.GLMeshItem(meshdata=gl.MeshData(vertexes=self.drop_dirs, faces=self.drop_faces), smooth=True,
                                  shader="shaded")
        self.bead_m = Arrow(v, WHITE, 3.0)
        for it in [self.drop0, self.drop1, self.bead] + self.drop_labels:
            v.addItem(it)

    def _build_curve(self):
        p = self.panels["curve"]; v = p.view
        self.crv_full = gl.GLLinePlotItem(pos=np.zeros((2, 3)), color=(0.7, 0.75, 0.85, 0.18), width=1.5, antialias=True)
        self.crv_grown = gl.GLLinePlotItem(pos=np.zeros((2, 3)), width=3.0, antialias=True)
        self.crv_start = gl.GLScatterPlotItem(pos=np.zeros((1, 3)), size=11, color=(0.25, 0.85, 0.35, 1), pxMode=True)
        self.crv_now = gl.GLScatterPlotItem(pos=np.zeros((1, 3)), size=10, color=WHITE, pxMode=True)
        for it in (self.crv_full, self.crv_grown, self.crv_start, self.crv_now):
            v.addItem(it)
        self.crv_tangent = Arrow(v, CYAN, 3.0, head=0.08)
        for k in range(3):                                                       # small x / y / z axes
            e = np.zeros(3); e[k] = 0.3
            c = FRAME_COLORS[k]
            v.addItem(gl.GLLinePlotItem(pos=np.array([np.zeros(3), e]), color=(c[0], c[1], c[2], 0.6), width=1.5))

    # ------------------------------------------------------------- data ---
    def ensemble(self):
        offs = np.array(parse_values(self.off_edit.text(), "ensemble offsets")) * 1e3
        b1 = np.array(parse_values(self.b1_edit.text(), "B1")) / 100.0
        if np.any(b1 <= 0):
            raise ValueError("B1 values must be positive")
        if offs.size * b1.size > 400:
            raise ValueError(f"{offs.size * b1.size} spins: at most 400")
        return offs, b1

    def pulse_changed(self, ps=None):
        """The main window's pulse changed: recompute when visible (else when shown)."""
        self._dirty = True
        if self.isVisible():
            self.recompute()

    def showEvent(self, event):
        super().showEvent(event)
        if self._dirty and HAVE_GL:
            self.recompute()

    def hideEvent(self, event):
        super().hideEvent(event)
        if HAVE_GL:
            self.pause()

    def _schedule(self, *args):
        self._recompute_timer.start()

    def recompute(self):
        if not HAVE_GL:
            return
        ps = self.pulse_provider() if self.pulse_provider else None
        if ps is None or ps.N == 0:
            self.error.setText("no pulse in the main window"); return
        try:
            offs, b1 = self.ensemble()
            frac = self.tau / self.dyn.N if self.dyn is not None and self.dyn.N else 0.0
            self.dyn = Dynamics(ps, offs, b1, [s.value() for s in self.m0_spins],
                                (self.sel_off.value() * 1e3, self.sel_b1.value() / 100))
        except ValueError as e:
            self.error.setText(str(e)); return
        self.error.setText("")
        self._dirty = False
        d = self.dyn
        self.tau = frac * d.N
        t = np.arange(d.N + 1) * ps.dt * 1e6
        self.amp_curve.setData(t, np.asarray(ps.amplitude, float))
        self.pulse_plot.setXRange(0, t[-1], padding=0.01)
        lo, hi = float(np.min(d.offsets[:-1])), float(np.max(d.offsets[:-1]))
        self._ens_colors = offset_colors(d.offsets[:-1], lo, hi)
        s_lo, s_hi = float(np.min(d.scales[:-1])), float(np.max(d.scales[:-1]))
        self._ens_sizes = (12 + 10 * (d.scales[:-1] - s_lo) / (s_hi - s_lo)) if s_hi > s_lo else np.full(d.n_ensemble, 15.0)
        self.colorbar.set_labels(f"{lo / 1e3:g} kHz", f"{hi / 1e3:g} kHz")
        b1txt = f"size: B1 {s_lo * 100:g} … {s_hi * 100:g} %" if s_hi > s_lo else f"B1 {s_lo * 100:g} %"
        self.panels["ensemble"].footer.setText(f"{d.n_ensemble} spins · colour: offset · {b1txt} · ringed: selected spin")
        c = d.curve
        self._crv_center = (c.min(axis=0) + c.max(axis=0)) / 2
        self._crv_scale = 2.0 / max(float(np.max(c.max(axis=0) - c.min(axis=0))), 1e-9)
        self.crv_full.setData(pos=(c - self._crv_center) * self._crv_scale)
        self.crv_start.setData(pos=((c[:1] - self._crv_center) * self._crv_scale))
        qf = d.q[-1, -1]
        for k, gh in enumerate(self.rot_ghosts):
            e = np.zeros(3); e[k] = 1.0
            gh.setData(pos=np.array([np.zeros(3), rotate(qf, e)]))
        self.draw()

    # --------------------------------------------------------- playback ---
    def toggle_play(self):
        if self.timer.isActive():
            self.pause()
        else:
            if self.dyn is None:
                self.recompute()
            if self.dyn is None:
                return
            if self.tau >= self.dyn.N - 1e-9:
                self.tau = 0.0
            import time
            self._last = time.perf_counter()
            self.timer.start(); self.play_button.setText("⏸ Pause")

    def pause(self):
        self.timer.stop(); self.play_button.setText("▶ Play")

    def _tick(self):
        import time
        now = time.perf_counter()
        dt = now - (self._last or now); self._last = now
        d = self.dyn
        if d is None:
            return self.pause()
        self.tau += d.N * dt / self.speed_spin.value()
        if self.tau >= d.N:
            if self.loop_check.isChecked():
                self.tau = 0.0
            else:
                self.tau = float(d.N); self.pause()
        self.draw()

    def set_time(self, tau):
        self.tau = float(np.clip(tau, 0, self.dyn.N if self.dyn else 0))
        self.draw()

    def _slider_moved(self, v):
        if self.dyn is not None and not getattr(self, "_syncing", False):
            self.set_time(v / 1000 * self.dyn.N)

    def _cursor_moved(self):
        if self.dyn is not None and not getattr(self, "_syncing", False):
            self.set_time(self.cursor.value() / (self.dyn.ps.dt * 1e6))

    # ------------------------------------------------------------- draw ---
    def _rot_mode_changed(self, *args):
        mode = self.rot_mode.currentIndex()
        for it in self.rot_frame_items:
            it.setVisible(mode == 0)
        for a in self.rot_arrows:
            a.setVisible(mode == 0)
        for it in [self.drop0, self.drop1] + self.drop_labels:
            it.setVisible(mode == 1)
        self.bead.setVisible(mode == 2); self.bead_m.setVisible(mode == 2)
        self.panels["rotation"].footer.setText(
            ("x / y / z rotated by R(t); faint: the final frame" if mode == 0 else
             "DROPS of U(t): top rank 0 (identity part), bottom rank 1 (rotation-axis part); colour = phase "
             "(red 0, yellow 90°, green 180°, blue 270°)" if mode == 1 else
             "Q-bead: sphere coloured by m · r (red: along m, green: opposite)"))
        self.draw()

    TEXT_INTERVAL = 0.2              # s: 2D widgets (slider, cursor, labels) during playback; each of their
    #                                  updates makes the window recomposite the GL panels (slow on macOS)

    def draw(self):
        d = self.dyn
        if d is None or not HAVE_GL:
            return
        tau = self.tau
        import time
        now = time.perf_counter()
        self._texts = not self.timer.isActive() or now - getattr(self, "_text_time", 0.0) >= self.TEXT_INTERVAL
        if self._texts:
            self._text_time = now
            self._syncing = True
            self.slider.setValue(int(round(1000 * tau / d.N)) if d.N else 0)
            self.cursor.setValue(tau * d.ps.dt * 1e6)
            self._syncing = False
            self.time_label.setText(f"t = {tau * d.ps.dt * 1e6:.2f} / {d.ps.duration * 1e6:.2f} µs")
        m = d.magnetization(tau)
        trail_len = int(round(self.trail_spin.value() / 100 * d.N))
        if self.panels["ensemble"].isVisible():
            self._draw_ensemble(d, tau, m, trail_len)
        if self.panels["field"].isVisible():
            self._draw_field(d, tau, m, trail_len)
        if self.panels["rotation"].isVisible():
            self._draw_rotation(d, tau)
        if self.panels["curve"].isVisible():
            self._draw_curve(d, tau)

    def _draw_ensemble(self, d, tau, m, trail_len):
        E = d.n_ensemble
        self.ens_dots.setData(pos=m[:E], color=self._ens_colors, size=self._ens_sizes)
        self.ens_sel.setData(pos=m[-1:])
        if trail_len > 0:
            tr = d.trail(tau, trail_len)[:E]                    # (E, L, 3)
            L = tr.shape[1]
            if L >= 2:
                seg = np.stack([tr[:, :-1], tr[:, 1:]], axis=2).reshape(-1, 3)
                fade = np.repeat(np.linspace(0.05, 0.9, L - 1), 2)
                col = np.repeat(self._ens_colors[:, None, :], 2 * (L - 1), axis=1).copy()
                col[:, :, 3] = fade[None]
                self.ens_trails.setData(pos=seg, color=col.reshape(-1, 4))
            self.ens_trails.setVisible(L >= 2)
        else:
            self.ens_trails.setVisible(False)
        if self.target_check.isChecked():
            tv = np.array([s.value() for s in self.target_spins])
            if np.linalg.norm(tv) > 0:
                tv /= np.linalg.norm(tv)
                self.ens_target.setData(pos=_circle(tv, tv * np.cos(0.12), np.sin(0.12), 48))
            self.ens_target.setVisible(np.linalg.norm(tv) > 0)
        else:
            self.ens_target.setVisible(False)

    def _draw_field(self, d, tau, m, trail_len):
        ms = m[-1]
        b = d.effective_field(tau)
        nb = np.linalg.norm(b)
        bh = b / nb if nb > 0 else np.array([0.0, 0.0, 1.0])
        self.fld_m.set(np.zeros(3), ms)
        self.fld_b.set(np.zeros(3), bh * 1.15)
        c = ms @ bh
        r = np.sqrt(max(1 - c * c, 0.0))
        self.fld_circle.setData(pos=_circle(bh, bh * c, r) if r > 1e-6 else np.zeros((2, 3)))
        tr = d.trail(tau, max(trail_len, 1))[-1]
        col = np.column_stack([np.tile(CYAN[:3], (len(tr), 1)), np.linspace(0.05, 0.9, len(tr))])
        self.fld_trail.setData(pos=tr, color=col)
        k = d.sub_pulse(tau)
        p = self.panels["field"]
        if self._texts:
            p.footer.setText(f"spin {d.offsets[-1] / 1e3:g} kHz, {d.scales[-1] * 100:g} % · sub-pulse {k + 1}: "
                         f"|B_eff| = {nb / 1e3:.2f} kHz (B1 {np.hypot(b[0], b[1]) / 1e3:.2f} kHz, "
                             f"offset {b[2] / 1e3:.2f} kHz) · cyan m, orange B_eff")

    def _draw_rotation(self, d, tau):
        q = d.quaternions(tau)[-1]
        mode = self.rot_mode.currentIndex()
        if mode == 0:
            k1 = int(np.floor(tau))
            path = d.q[-1, :k1 + 1]
            for k, (a, tip) in enumerate(zip(self.rot_arrows, self.rot_tips)):
                e = np.zeros(3); e[k] = 1.0
                a.set(np.zeros(3), rotate(q, e))
                pts = np.vstack([rotate(path, e), rotate(q, e)])
                tip.setData(pos=pts)
        elif mode == 1:
            f0, f1 = propagator_droplets(q, self.drop_dirs)
            for item, f, z in ((self.drop0, f0, 0.8), (self.drop1, f1, -0.7)):     # stacked: fit tall panels
                verts = droplet(f, self.drop_dirs) * 0.6 + np.array([0, 0, z])
                item.setMeshData(meshdata=mesh(verts, self.drop_faces, phase_colors(f)))
        else:
            ms = rotate(q, d.m0)
            f = magnetization_droplet(ms, self.drop_dirs).real
            col = np.empty((len(f), 4))
            w = np.abs(f)[:, None]
            base = np.array([0.85, 0.87, 0.90])
            col[:, :3] = np.where(f[:, None] >= 0, base + (np.array([0.92, 0.18, 0.18]) - base) * w,
                                  base + (np.array([0.20, 0.78, 0.28]) - base) * w)
            col[:, 3] = 1.0
            self.bead.setMeshData(meshdata=mesh(self.drop_dirs * 0.85, self.drop_faces, col))
            self.bead_m.set(np.zeros(3), ms * 1.3)

    def _draw_curve(self, d, tau):
        pts = (d.curve_until(tau) - self._crv_center) * self._crv_scale
        self.crv_grown.setData(pos=pts, color=time_colors(len(pts))) if len(pts) >= 2 else None
        self.crv_now.setData(pos=pts[-1:])
        qs = d.quaternions(tau)[-1].copy(); qs[1:] *= -1
        t = rotate(qs, [0.0, 0.0, 1.0])
        self.crv_tangent.set(pts[-1], t * 0.3)
        if self._texts:
            self.panels["curve"].footer.setText(f"spin {d.offsets[-1] / 1e3:g} kHz, {d.scales[-1] * 100:g} % · "
                                                "green → red: time · cyan: tangent (toggling-frame z)")

    # ----------------------------------------------------------- export ---
    def frame_image(self):
        """The visible panels side by side as one QImage (None if none is visible)."""
        if not self.mv.visible_keys() or not self.mv.gl_ok():
            return None
        return self.mv.grabFramebuffer()

    def export_gif(self, path, fps=None, max_frames=600):
        """One pass of the pulse in the playback time as a GIF at ``fps`` frames per second."""
        from PIL import Image
        from PySide6.QtWidgets import QApplication
        if self.dyn is None:
            raise ValueError("nothing to export")
        fps = fps or self.export_fps.value()
        n = int(min(max_frames, max(2, round(self.speed_spin.value() * fps))))
        was = self.timer.isActive(); self.pause()
        frames = []
        for i in range(n):
            self.set_time(i / (n - 1) * self.dyn.N)
            QApplication.processEvents()
            img = self.frame_image()
            if img is None:
                break
            img = img.convertToFormat(QImage.Format.Format_RGBA8888)
            buf = img.constBits()
            arr = np.frombuffer(buf, np.uint8, count=img.sizeInBytes()).reshape(img.height(), img.bytesPerLine())
            arr = arr[:, :img.width() * 4].reshape(img.height(), img.width(), 4)
            frames.append(Image.fromarray(arr[:, :, :3].copy()).quantize(colors=255))
        if not frames:
            raise ValueError("no panel visible")
        frames[0].save(path, save_all=True, append_images=frames[1:], duration=int(round(1000 / fps)), loop=0)
        if was:
            self.toggle_play()
        return len(frames)

    def export_gif_dialog(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export GIF", "pulse_animation.gif", "GIF (*.gif)")
        if not path:
            return
        try:
            n = self.export_gif(path)
            self.error.setStyleSheet("color: #4ade80;")
            self.error.setText(f"saved {n} frames: {path}")
        except Exception as e:
            self.error.setStyleSheet("color: #f87171;")
            self.error.setText(f"export failed: {e}")

    def shutdown(self):
        if HAVE_GL:
            self.pause()
