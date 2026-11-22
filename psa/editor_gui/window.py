"""Pulse-sequence editor window: tabs (Text, Waveform, ...) and the shared summary panel."""
from __future__ import annotations

import os

from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QFrame, QHBoxLayout, QMainWindow, QScrollArea, QSplitter, QTabWidget, QWidget

from psa.editor_gui.panel import SummaryPanel
from psa.editor_gui.tab_text import TextTab
from psa.editor_gui.tab_waveform import WaveformTab
from psa.editor_gui.tab_families import FamiliesTab


def scrollable(widget):
    """Wrap ``widget`` in a vertically scrollable area (content keeps its minimum height)."""
    area = QScrollArea()
    area.setWidgetResizable(True)
    area.setFrameShape(QFrame.Shape.NoFrame)
    area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    area.setWidget(widget)
    # only vertical scrolling: reserve the content's minimum width (+ room for the scroll bar)
    area.setMinimumWidth(widget.minimumSizeHint().width() + area.verticalScrollBar().sizeHint().width() + 4)
    return area


class PulseEditorWindow(QMainWindow):
    """``main_window`` (optional) is the PySide6 analyzer; "Save + plot" loads the saved file there."""

    def __init__(self, main_window=None, dt=0.5e-6, umax=1e4, save_path=None, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("Pulse Sequence Editor")
        save_path = save_path or os.path.join(os.path.expanduser("~"), "Edited_Pulse_Sequence.csv")

        self.panel = SummaryPanel(dt, umax, save_path)
        self.tabs = QTabWidget()
        self.text_tab = TextTab(dt, umax)
        self.waveform_tab = WaveformTab(dt, umax)
        self.families_tab = FamiliesTab(dt, umax)
        self._tab_list = []                      # inner tab widgets, in tab order
        for tab, title in ((self.text_tab, "Text"), (self.waveform_tab, "Waveform"),
                           (self.families_tab, "Curve families")):
            self._tab_list.append(tab)
            self.tabs.addTab(scrollable(tab), title)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self.tabs)
        splitter.addWidget(scrollable(self.panel))
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)
        central = QWidget()
        QHBoxLayout(central).addWidget(splitter)
        self.setCentralWidget(central)
        self._fit_to_screen()

        for tab in self.all_tabs():
            tab.pulseChanged.connect(lambda ps, t=tab: self._pulse_changed(t, ps))
        self.waveform_tab.requestTextPulse.connect(self._waveform_from_text)
        self.families_tab.editInWaveform.connect(self._edit_in_waveform)
        self.tabs.currentChanged.connect(self._tab_switched)
        self.panel.metadataChanged.connect(self._metadata_changed)
        self.panel.sendRequested.connect(self._send_to_main)
        self.text_tab.parse()
        self._tab_switched(self.tabs.currentIndex())

    def _fit_to_screen(self):
        screen = QGuiApplication.primaryScreen()
        if screen is None:
            self.resize(1500, 900)
            return
        avail = screen.availableGeometry()
        self.resize(min(1500, int(avail.width() * 0.95)), min(950, int(avail.height() * 0.92)))

    def all_tabs(self):
        return list(self._tab_list)

    def current_tab(self):
        i = self.tabs.currentIndex()
        return self._tab_list[i] if 0 <= i < len(self._tab_list) else None

    def show_tab(self, tab):
        self.tabs.setCurrentIndex(self._tab_list.index(tab))

    # --------------------------------------------------------------- events ---
    def _pulse_changed(self, tab, ps):
        if tab is not self.current_tab():
            return
        if ps is None:                            # e.g. live calculation switched off
            self.panel.set_pulse(None)
            return
        if abs(ps.dt - self.panel.dt) > 1e-15 or abs(ps.umax - self.panel.umax) > 1e-9:   # e.g. a loaded design
            for spin, v in ((self.panel.dt_spin, ps.dt * 1e6), (self.panel.umax_spin, ps.umax / 1e3)):
                spin.blockSignals(True); spin.setValue(v); spin.blockSignals(False)
            for t in self.all_tabs():
                if t is not tab:
                    t.set_metadata(ps.dt, ps.umax)
        self.panel.set_pulse(ps)

    def _tab_switched(self, index):
        tab = self.current_tab()
        if hasattr(tab, "current_pulse"):
            self.panel.set_pulse(tab.current_pulse())

    def _metadata_changed(self, dt, umax):
        for tab in self.all_tabs():
            tab.set_metadata(dt, umax)
        self._tab_switched(self.tabs.currentIndex())

    def _waveform_from_text(self):
        ps = self.text_tab.current_pulse()
        if ps is not None:
            self.waveform_tab.load_pulse(ps)

    def _edit_in_waveform(self, ps):
        self.waveform_tab.load_pulse(ps)
        self.show_tab(self.waveform_tab)

    def _send_to_main(self, path):
        mw = self.main_window
        if mw is None:
            return
        ps = self.panel.pulse
        mw.pulse_sequence_edit.setText(path)
        mw.time_per_pulse_spin.setValue(ps.dt * 1e6)
        mw.pulse_amplitude_spin.setValue(ps.umax / 1e3)
        mw.update_analysis()
        if hasattr(mw, "append_info"):
            mw.append_info(f"Pulse from editor loaded: {path}")
