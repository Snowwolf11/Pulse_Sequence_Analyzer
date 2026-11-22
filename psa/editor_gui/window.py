"""Pulse-sequence editor window: tabs (Text, Waveform, ...) and the shared summary panel."""
from __future__ import annotations

import os

from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (QDialog, QFrame, QHBoxLayout, QMainWindow, QPushButton, QScrollArea, QSplitter,
                               QTabWidget, QTextEdit, QVBoxLayout, QWidget)

from psa.editor_gui.panel import SummaryPanel
from psa.editor_gui.widgets import lock_heights
from psa.editor_gui.tab_text import TextTab
from psa.editor_gui.tab_waveform import WaveformTab
from psa.editor_gui.tab_families import FamiliesTab
from psa.editor_gui.tab_builder import BuilderTab
from psa.editor_gui.tab_import import ImportTab
from psa.editor_gui.tab_grape import GrapeTab


class FitScrollArea(QScrollArea):
    """Scroll area that tells its widget the visible height (``fit_to_height``), so that tabs can keep
    their main plots as large as the visible area and put extra content (e.g. the solver) below."""

    def resizeEvent(self, event):
        super().resizeEvent(event)
        w = self.widget()
        if hasattr(w, "fit_to_height"):
            w.fit_to_height(self.viewport().height())


def scrollable(widget):
    """Wrap ``widget`` in a vertically scrollable area (content keeps its minimum height)."""
    area = FitScrollArea()
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

        from psa.design.optimize import SolverSettings
        self.solver_settings = getattr(main_window, "solver_settings", None) or SolverSettings()
        self.panel = SummaryPanel(dt, umax, save_path)
        self.tabs = QTabWidget()
        self.text_tab = TextTab(dt, umax)
        self.waveform_tab = WaveformTab(dt, umax)
        self.families_tab = FamiliesTab(dt, umax)
        self.builder_tab = BuilderTab(dt, umax)
        self.import_tab = ImportTab(dt, umax)
        self.grape_tab = GrapeTab(dt, umax)
        self._tab_list = []                      # inner tab widgets, in tab order
        for tab, title in ((self.text_tab, "Text"), (self.waveform_tab, "Waveform"),
                           (self.families_tab, "Curve families"), (self.builder_tab, "Curve builder"),
                           (self.import_tab, "Import curve"), (self.grape_tab, "Optimize (GRAPE)")):
            self._tab_list.append(tab)
            self.tabs.addTab(scrollable(tab), title)

        self.help_button = QPushButton("Help")
        self.help_button.clicked.connect(self.show_help)
        self.tabs.setCornerWidget(self.help_button, Qt.Corner.TopRightCorner)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self.tabs)
        splitter.addWidget(scrollable(self.panel))
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)
        central = QWidget()
        QHBoxLayout(central).addWidget(splitter)
        self.setCentralWidget(central)
        lock_heights(self)
        self._fit_to_screen()

        for tab in self.all_tabs():
            tab.pulseChanged.connect(lambda ps, t=tab: self._pulse_changed(t, ps))
        for tab in self.all_tabs():
            if hasattr(tab, "solver_box"):
                tab.solver_box.settings_provider = lambda: self.solver_settings
                tab.solver_box.settingsRequested.connect(self.open_solver_settings)
        self.grape_tab.settings_provider = lambda: self.solver_settings
        self.grape_tab.settingsRequested.connect(self.open_solver_settings)
        self.grape_tab.requestStart.connect(self._grape_start)
        self.grape_tab.editInWaveform.connect(self._edit_in_waveform)
        self.waveform_tab.requestTextPulse.connect(self._waveform_from_text)
        self.families_tab.editInWaveform.connect(self._edit_in_waveform)
        self.builder_tab.editInWaveform.connect(self._edit_in_waveform)
        self.import_tab.editInWaveform.connect(self._edit_in_waveform)
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

    def open_solver_settings(self):
        """Solver tab of the analyzer's Settings dialog (or a stand-alone dialog without analyzer)."""
        mw = self.main_window
        if mw is not None and hasattr(mw, "solver_settings"):
            from psa.GUI_pyside6 import SettingsDialog
            self._settings_dialog = SettingsDialog(mw, initial_tab="Solver")
        else:
            from psa.editor_gui.solver_settings import SolverSettingsDialog
            self._settings_dialog = SolverSettingsDialog(self.solver_settings, self)
        self._settings_dialog.exec()

    def show_help(self):
        from psa.documentation import EDITOR_DOCUMENTATION_STR
        if getattr(self, "_help", None) is None:
            self._help = QDialog(self)
            self._help.setWindowTitle("Pulse Sequence Editor – Help")
            self._help.resize(820, 760)
            lay = QVBoxLayout(self._help)
            text = QTextEdit(); text.setReadOnly(True); text.setPlainText(EDITOR_DOCUMENTATION_STR)
            lay.addWidget(text)
            close = QPushButton("Close"); close.clicked.connect(self._help.close)
            lay.addWidget(close)
        self._help.show(); self._help.raise_()

    def closeEvent(self, event):
        for tab in self.all_tabs():
            if hasattr(tab, "solver_box"):
                tab.solver_box.stop(wait=True)
        self.grape_tab.stop(wait=True)
        super().closeEvent(event)

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

    def _grape_start(self, source):
        tabs = {"Text tab": self.text_tab, "Waveform tab": self.waveform_tab, "Curve families tab": self.families_tab,
                "Curve builder tab": self.builder_tab, "Import curve tab": self.import_tab}
        ps = tabs[source].current_pulse() if source in tabs else None
        if ps is None:
            self.grape_tab.start_failed(f"the {source} has no valid pulse (or its live calculation is off)")
        else:
            self.grape_tab.set_start(ps)

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
