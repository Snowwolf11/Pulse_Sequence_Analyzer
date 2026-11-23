"""Small dialog to choose the goal of GRAPE points: the default (the GRAPE tab's gate goal), a UR rotation or a
PP transfer, with presets for selective pulses."""
from __future__ import annotations

from PySide6.QtWidgets import (QButtonGroup, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
                               QHBoxLayout, QLabel, QRadioButton, QVBoxLayout, QWidget)

from psa.editor_gui.widgets import BAD, MUTED_FG, lock_heights

PRESETS = (
    ("presets …", None),
    ("do nothing (UR 0°)", {"gate": "UR", "axis": [0, 0, 1], "angle": 0.0}),
    ("excite z → x (PP)", {"gate": "PP", "m0": [0, 0, 1], "mt": [1, 0, 0]}),
    ("excite z → y (PP)", {"gate": "PP", "m0": [0, 0, 1], "mt": [0, 1, 0]}),
    ("invert z → −z (PP)", {"gate": "PP", "m0": [0, 0, 1], "mt": [0, 0, -1]}),
    ("keep z (PP z → z)", {"gate": "PP", "m0": [0, 0, 1], "mt": [0, 0, 1]}),
    ("UR 90° about x", {"gate": "UR", "axis": [1, 0, 0], "angle": 90.0}),
    ("UR 180° about x", {"gate": "UR", "axis": [1, 0, 0], "angle": 180.0}),
)


def _vector(values):
    row = QHBoxLayout(); spins = []
    for name, v in zip("xyz", values):
        sp = QDoubleSpinBox(); sp.setRange(-1, 1); sp.setDecimals(3); sp.setSingleStep(0.1); sp.setValue(v)
        sp.setPrefix(name + " ")
        row.addWidget(sp); spins.append(sp)
    return row, spins


class GoalDialog(QDialog):
    """Modal: ``GoalDialog.ask(parent, goal)`` -> (accepted, goal or None)."""

    def __init__(self, goal=None, count=1, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Goal of the selected points")
        lay = QVBoxLayout(self)
        intro = QLabel(f"{count} point{'s' * (count != 1)} selected. Default: the gate goal of the GRAPE tab. "
                       "Own goals make selective pulses possible (e.g. excite inside a band, do nothing outside).")
        intro.setWordWrap(True); intro.setStyleSheet(f"color: {MUTED_FG};")
        lay.addWidget(intro)
        self.default_radio = QRadioButton("default: the gate goal of the GRAPE tab")
        self.ur_radio = QRadioButton("UR rotation (the whole rotation)")
        self.pp_radio = QRadioButton("PP transfer (m0 → target)")
        g = QButtonGroup(self)
        for r in (self.default_radio, self.ur_radio, self.pp_radio):
            g.addButton(r)
        self.preset_box = QComboBox(); self.preset_box.addItems([p[0] for p in PRESETS])
        self.preset_box.activated.connect(self._preset)
        lay.addWidget(self.preset_box)
        lay.addWidget(self.default_radio)
        lay.addWidget(self.ur_radio)
        self.ur_widget = QWidget(); uf = QFormLayout(self.ur_widget); uf.setContentsMargins(18, 0, 0, 0)
        self.angle_spin = QDoubleSpinBox(); self.angle_spin.setRange(-720, 720); self.angle_spin.setDecimals(3)
        self.angle_spin.setSuffix(" °"); self.angle_spin.setValue(90.0)
        row, self.axis_spins = _vector((1.0, 0.0, 0.0))
        uf.addRow("angle:", self.angle_spin); uf.addRow("axis:", row)
        lay.addWidget(self.ur_widget)
        lay.addWidget(self.pp_radio)
        self.pp_widget = QWidget(); pf = QFormLayout(self.pp_widget); pf.setContentsMargins(18, 0, 0, 0)
        row, self.m0_spins = _vector((0.0, 0.0, 1.0)); pf.addRow("from m0:", row)
        row, self.mt_spins = _vector((1.0, 0.0, 0.0)); pf.addRow("to:", row)
        lay.addWidget(self.pp_widget)
        self.error = QLabel(""); self.error.setStyleSheet(f"color: {BAD};")
        lay.addWidget(self.error)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self._ok); buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)
        for r in (self.default_radio, self.ur_radio, self.pp_radio):
            r.toggled.connect(self._kind_changed)
        self.set_goal(goal)
        lock_heights(self)

    def set_goal(self, goal):
        if goal is None:
            self.default_radio.setChecked(True)
        elif goal["gate"] == "UR":
            self.ur_radio.setChecked(True)
            self.angle_spin.setValue(float(goal["angle"]))
            for sp, v in zip(self.axis_spins, goal["axis"]):
                sp.setValue(float(v))
        else:
            self.pp_radio.setChecked(True)
            for sp, v in zip(self.m0_spins + self.mt_spins, list(goal["m0"]) + list(goal["mt"])):
                sp.setValue(float(v))
        self._kind_changed()

    def _preset(self, i):
        if PRESETS[i][1] is not None:
            self.set_goal(PRESETS[i][1])
        self.preset_box.setCurrentIndex(0)

    def _kind_changed(self, *args):
        self.ur_widget.setEnabled(self.ur_radio.isChecked())
        self.pp_widget.setEnabled(self.pp_radio.isChecked())

    def goal(self):
        """The chosen goal (None = default); ValueError if a vector is zero."""
        if self.default_radio.isChecked():
            return None
        if self.ur_radio.isChecked():
            goal = {"gate": "UR", "axis": [s.value() for s in self.axis_spins], "angle": self.angle_spin.value()}
        else:
            goal = {"gate": "PP", "m0": [s.value() for s in self.m0_spins], "mt": [s.value() for s in self.mt_spins]}
        from psa.design.grape import check_goal
        check_goal(goal)
        return goal

    def _ok(self):
        try:
            self.goal()
        except ValueError as e:
            self.error.setText(str(e))
            return
        self.accept()

    @classmethod
    def ask(cls, parent, goal=None, count=1):
        dlg = cls(goal, count, parent)
        ok = dlg.exec() == QDialog.DialogCode.Accepted
        return ok, (dlg.goal() if ok else None)
