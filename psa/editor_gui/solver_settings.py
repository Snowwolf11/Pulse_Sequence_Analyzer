"""Solver settings widget (used in the analyzer's Settings dialog and from the editor)."""
from __future__ import annotations

from PySide6.QtWidgets import (QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QGroupBox, QLabel, QLineEdit,
                               QPushButton, QSpinBox, QVBoxLayout, QWidget)

from psa.design.optimize import SolverSettings
from psa.editor_gui.widgets import BAD, MUTED_FG, lock_heights

SCI_FIELDS = (("xtol", "step tolerance (xtol)", "stop when the step is this small (relative to the variables); 0 = off"),
              ("ftol", "improvement tolerance (ftol)", "stop when the residual improves by less than this (relative); 0 = off"),
              ("gtol", "gradient tolerance (gtol)", "stop when the gradient is this small; 0 = off"),
              ("diff_step", "finite-difference step", "relative step for the numerical derivatives (0 = automatic, ≈1.5e-8)"),
              ("solved_tolerance", "'solved' below", "a robustness goal (|A_k|/T^k, |B1|/∫Ω) counts as reached below this"),
              ("gate_tolerance_deg", "'solved' gate error below [°]", "the target gate counts as reached below this error "
               "(Waveform, Curve families, Curve builder)"),
              ("gate_infidelity_tolerance", "'solved' gate infidelity below (GRAPE)",
               "Optimize (GRAPE): the gate goal counts as reached when the infidelity is below this at every point "
               "(1e-12 ≈ 1.4e-4° rotation error)"))
WEIGHTS = (("weight_order1", "1st order"), ("weight_order2", "2nd order"), ("weight_order3", "3rd order"),
           ("weight_gate", "target gate"), ("weight_b1", "B1 scaling"))


class SolverSettingsWidget(QWidget):
    """Edits a SolverSettings object; ``apply()`` writes the fields back (returns an error text or "")."""

    def __init__(self, settings: SolverSettings, parent=None):
        super().__init__(parent)
        self.settings = settings
        layout = QVBoxLayout(self)
        intro = QLabel("These settings are used by the solver in the pulse sequence editor (Waveform, "
                       "Curve families, Curve builder and Optimize (GRAPE) tabs). Goals and free variables are chosen in the tabs.")
        intro.setWordWrap(True); intro.setStyleSheet(f"color: {MUTED_FG};")
        layout.addWidget(intro)

        stop = QGroupBox("Stopping")
        sf = QFormLayout(stop)
        self.solved_check = QCheckBox("stop as soon as all goals are reached")
        self.solved_check.setToolTip("on: stop when every goal is below its 'solved' threshold (below); "
                                     "off: keep refining until the solver converges")
        sf.addRow(self.solved_check)
        self.max_spin = QSpinBox(); self.max_spin.setRange(10, 100_000_000)
        self.max_spin.setToolTip("maximum number of evaluations in total (including those for the numerical "
                                 "derivatives; the counter in the solver box)")
        sf.addRow("max. evaluations:", self.max_spin)
        self.time_spin = QDoubleSpinBox(); self.time_spin.setRange(0, 1e7); self.time_spin.setDecimals(1)
        self.time_spin.setSuffix(" s"); self.time_spin.setToolTip("stop after this time and keep the best result (0 = no limit)")
        sf.addRow("time limit (0 = none):", self.time_spin)
        self.sci = {}
        for name, label, tip in SCI_FIELDS:
            e = QLineEdit(); e.setToolTip(tip)
            sf.addRow(label + ":", e)
            self.sci[name] = e
        layout.addWidget(stop)

        wbox = QGroupBox("Goal weights (relative importance in the least-squares residual)")
        wf = QFormLayout(wbox)
        self.weights = {}
        for name, label in WEIGHTS:
            sp = QDoubleSpinBox(); sp.setRange(0.0, 1e6); sp.setDecimals(4)
            wf.addRow(label + ":", sp)
            self.weights[name] = sp
        layout.addWidget(wbox)

        gbox = QGroupBox("GRAPE point editor: lock to grid")
        gf = QFormLayout(gbox)
        self.grid_off_spin = QDoubleSpinBox(); self.grid_off_spin.setRange(0, 1e4); self.grid_off_spin.setDecimals(3)
        self.grid_off_spin.setSuffix(" kHz"); self.grid_off_spin.setSpecialValueText("off")
        self.grid_b1_spin = QDoubleSpinBox(); self.grid_b1_spin.setRange(0, 100); self.grid_b1_spin.setDecimals(3)
        self.grid_b1_spin.setSuffix(" %"); self.grid_b1_spin.setSpecialValueText("off")
        for sp in (self.grid_off_spin, self.grid_b1_spin):
            sp.setToolTip("with 'lock to grid' on, points snap to multiples of this step; off = 0")
        gf.addRow("offset step:", self.grid_off_spin)
        gf.addRow("B1 step:", self.grid_b1_spin)
        layout.addWidget(gbox)

        self.error = QLabel(""); self.error.setStyleSheet(f"color: {BAD};"); self.error.setWordWrap(True)
        layout.addWidget(self.error)
        reset = QPushButton("Reset solver settings to defaults"); reset.clicked.connect(lambda: self.load(SolverSettings()))
        layout.addWidget(reset)
        layout.addStretch()
        self.load(settings)
        lock_heights(self)

    def load(self, s: SolverSettings):
        self.max_spin.setValue(int(s.max_nfev))
        self.time_spin.setValue(float(s.time_limit_s))
        self.solved_check.setChecked(bool(s.stop_when_solved))
        for name, e in self.sci.items():
            e.setText(f"{getattr(s, name):.3g}")
        for name, sp in self.weights.items():
            sp.setValue(getattr(s, name))
        self.grid_off_spin.setValue(float(s.point_grid_offset_khz)); self.grid_b1_spin.setValue(float(s.point_grid_b1_percent))
        self.error.setText("")

    def apply(self) -> str:
        values = {"max_nfev": self.max_spin.value(), "time_limit_s": self.time_spin.value(),
                  "stop_when_solved": self.solved_check.isChecked()}
        for name, e in self.sci.items():
            try:
                v = float(e.text().replace(",", "."))
            except ValueError:
                return self._fail(f"'{e.text()}' is not a number ({name})")
            if v < 0:
                return self._fail(f"{name} must not be negative")
            values[name] = v
        for name, sp in self.weights.items():
            values[name] = sp.value()
        values["point_grid_offset_khz"] = self.grid_off_spin.value()
        values["point_grid_b1_percent"] = self.grid_b1_spin.value()
        if not any(values[n] > 0 for n, _ in WEIGHTS):
            return self._fail("at least one weight must be positive")
        for k, v in values.items():
            setattr(self.settings, k, v)
        self.error.setText("")
        return ""

    def _fail(self, text):
        self.error.setText(text)
        return text


class SolverSettingsDialog(QDialog):
    """Stand-alone dialog (editor without analyzer window)."""

    def __init__(self, settings: SolverSettings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Solver settings")
        layout = QVBoxLayout(self)
        self.widget = SolverSettingsWidget(settings)
        layout.addWidget(self.widget)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self._ok); buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _ok(self):
        if not self.widget.apply():
            self.accept()
