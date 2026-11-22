"""Text tab: the existing R(...) / PR(...) pulse syntax, parsed live."""
from __future__ import annotations

from PySide6.QtCore import QTimer, Signal
from PySide6.QtWidgets import QLabel, QTextEdit, QVBoxLayout, QWidget

from psa.design.pulse import PulseSequence
from psa.editor_gui.widgets import BAD, MUTED_FG
from psa.pulse_string_to_array import parse_pulse_sequence

HELP = ("R(B1 amplitude [%], phase [°], duration [µs]);   PR(B1 amplitude [%], start phase [°], end phase [°], "
        "duration [µs]);   durations must be multiples of the sub-pulse duration")
EXAMPLE = "R(100, 0, 25);\nPR(100, 0, 360, 100);\nR(100, 180, 25);"


class TextTab(QWidget):
    pulseChanged = Signal(object)

    def __init__(self, dt=0.5e-6, umax=1e4, parent=None):
        super().__init__(parent)
        self.dt, self.umax = dt, umax
        self.pulse = None
        layout = QVBoxLayout(self)
        hint = QLabel(HELP); hint.setWordWrap(True); hint.setStyleSheet(f"color: {MUTED_FG};")
        layout.addWidget(hint)
        self.editor = QTextEdit()
        self.editor.setPlainText(EXAMPLE)
        self.editor.setMinimumHeight(250)
        layout.addWidget(self.editor)
        self.error = QLabel(""); self.error.setWordWrap(True); self.error.setStyleSheet(f"color: {BAD};")
        layout.addWidget(self.error)
        self._timer = QTimer(self); self._timer.setSingleShot(True); self._timer.setInterval(300)
        self._timer.timeout.connect(self.parse)
        self.editor.textChanged.connect(self._timer.start)

    def set_metadata(self, dt, umax):
        self.dt, self.umax = dt, umax
        self.parse()

    def set_text(self, text):
        self.editor.setPlainText(text)
        self.parse()

    def parse(self):
        try:
            PS = parse_pulse_sequence(self.editor.toPlainText(), subpulse_dt_us=self.dt * 1e6)
        except Exception as e:
            self.error.setText(str(e))
            return None
        if len(PS) == 0:
            self.error.setText("no R(...) or PR(...) blocks found")
            return None
        self.error.setText("")
        self.pulse = PulseSequence.from_array(PS, self.dt, self.umax)
        self.pulseChanged.emit(self.pulse)
        return self.pulse

    def current_pulse(self):
        return self.pulse if self.pulse is not None else self.parse()
