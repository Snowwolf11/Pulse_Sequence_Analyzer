"""Small window of a pulse library (a quick-mode series of durations): worst / mean infidelity vs. duration
and a table of the pulses; a click on a row loads that pulse. Can also open an existing library folder."""
from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QDialog, QFileDialog, QHBoxLayout, QLabel, QPushButton,
                               QTableWidget, QTableWidgetItem, QVBoxLayout)

from psa.design.library import read_index
from psa.design.pulse import PulseSequence
from psa.editor_gui.widgets import ACCENT, ACCENT_2, BAD, MUTED_FG, make_plot

COLUMNS = ("duration [µs]", "worst 1 − F", "mean 1 − F", "below 1e-4", "amp. rate [%/µs]", "phase rate [kHz]", "note")


class LibraryWindow(QDialog):
    pulseChosen = Signal(object)             # PulseSequence of the clicked row

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Pulse library")
        self.setWindowFlag(Qt.WindowType.Window, True)
        self.resize(760, 620)
        self.folder, self.pulses, self.total = None, [], 0
        lay = QVBoxLayout(self)
        top = QHBoxLayout()
        self.info = QLabel("no library yet"); self.info.setWordWrap(True); self.info.setStyleSheet(f"color: {MUTED_FG};")
        self.open_button = QPushButton("Open library…"); self.open_button.clicked.connect(self.choose_folder)
        self.open_button.setToolTip("show an existing library folder (with library.json)")
        top.addWidget(self.info, 1); top.addWidget(self.open_button)
        lay.addLayout(top)
        self.plot = make_plot("log10(1 − F)", "duration [µs]")
        self.plot.setMinimumHeight(260)
        self.plot.addLegend(offset=(-10, 10))
        self.worst_curve = self.plot.plot(pen=pg.mkPen(ACCENT, width=2), symbol="o", symbolSize=7,
                                          symbolBrush=ACCENT, name="worst")
        self.mean_curve = self.plot.plot(pen=pg.mkPen(ACCENT_2, width=2), symbol="o", symbolSize=7,
                                         symbolBrush=ACCENT_2, name="mean")
        lay.addWidget(self.plot, 1)
        self.table = QTableWidget(0, len(COLUMNS)); self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.cellClicked.connect(self._row_clicked)
        self.table.setToolTip("click a row to load that pulse")
        lay.addWidget(self.table, 1)
        self.error = QLabel(""); self.error.setStyleSheet(f"color: {BAD};")
        lay.addWidget(self.error)

    def set_library(self, folder, pulses, total=None):
        """Show a library: its folder and the entries of library.json (total = planned number of pulses)."""
        self.folder, self.pulses = folder, list(pulses)
        self.total = total if total is not None else len(self.pulses)
        self.error.setText("")
        self._refresh()

    def _refresh(self):
        ps = self.pulses
        name = os.path.basename(self.folder) if self.folder else "—"
        self.info.setText(f"{name}: {len(ps)} of {self.total} pulses" + (f"\n{self.folder}" if self.folder else ""))
        d = np.array([p["duration_us"] for p in ps], float)
        order = np.argsort(d) if len(d) else []
        for curve, key in ((self.worst_curve, "worst"), (self.mean_curve, "mean")):
            y = np.log10(np.maximum([ps[i][key] for i in order], 1e-16)) if len(d) else []
            curve.setData(d[order] if len(d) else [], y)
        self.table.setRowCount(len(ps))
        for r, i in enumerate(order):
            p = ps[i]
            note = "stopped" if p.get("stopped") else ("" if p.get("within_limits", True) else "rate limits broken")
            vals = (f"{p['duration_us']:g}", f"{p['worst']:.2e}", f"{p['mean']:.2e}", f"{100 * p['below_1e-4']:.0f} %",
                    f"{p['amp_rate']:.1f}", f"{p['phase_rate']:.1f}", note)
            for c, v in enumerate(vals):
                item = QTableWidgetItem(v)
                item.setData(Qt.ItemDataRole.UserRole, int(i))
                self.table.setItem(r, c, item)
        self.table.resizeColumnsToContents()

    def pulse(self, i):
        """PulseSequence of entry i (read from the library folder)."""
        return PulseSequence.load(os.path.join(self.folder, self.pulses[i]["file"]))

    def _row_clicked(self, row, col):
        item = self.table.item(row, 0)
        if item is None:
            return
        try:
            self.pulseChosen.emit(self.pulse(item.data(Qt.ItemDataRole.UserRole)))
        except (OSError, ValueError) as e:
            self.error.setText(f"cannot read the pulse: {e}")

    def choose_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Open a pulse library", self.folder or "")
        if folder:
            self.open_folder(folder)

    def open_folder(self, folder):
        try:
            d = read_index(folder)
        except (OSError, ValueError) as e:
            self.error.setText(f"not a pulse library: {e}")
            return False
        self.set_library(folder, d["pulses"], d.get("settings", {}).get("n_pulses"))
        return True
