"""Solver box shared by the waveform, curve-family and builder tabs.

Goals: detuning orders (1st/2nd/3rd), optional target gate, optional first-order B1 robustness.
The solver runs in a separate process (psa.design.optimize.run_job); the GUI stays responsive and
Stop returns the best variables found so far.
"""
from __future__ import annotations

import multiprocessing as mp
import queue as queue_module
import time

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (QCheckBox, QDoubleSpinBox, QGroupBox, QHBoxLayout, QLabel, QListWidget,
                               QListWidgetItem, QPushButton, QVBoxLayout)

from psa.design.optimize import Goals, run_job
from psa.editor_gui.widgets import BAD, GOOD, MUTED_FG, lock_heights


class SolverBox(QGroupBox):
    running = Signal(bool)
    solved = Signal(object, object)          # solution dict, info dict
    settingsRequested = Signal()

    def __init__(self, title="Solve", show_variables=True, variables_hint="", parent=None):
        super().__init__(title, parent)
        self.job_provider = None              # callable(goals, free_keys) -> job dict
        self.settings_provider = None         # callable() -> SolverSettings (None: defaults)
        self._process = None
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("cancel detuning error:"))
        orow = QHBoxLayout()
        self.order_checks = {k: QCheckBox(t) for k, t in ((1, "1st (closed)"), (2, "2nd (areas)"), (3, "3rd"))}
        self.order_checks[1].setChecked(True)
        for k in (1, 2, 3):
            orow.addWidget(self.order_checks[k])
        orow.addStretch()
        layout.addLayout(orow)
        self.b1_check = QCheckBox("B1 scaling robust (1st order)")
        self.b1_check.setToolTip("first-order robustness against amplitude (B1) scaling errors")
        layout.addWidget(self.b1_check)
        trow = QHBoxLayout()
        self.target_check = QCheckBox("target gate:")
        self.target_check.setToolTip("optional: the pulse must also implement this rotation (angle about axis)")
        self.angle_spin = QDoubleSpinBox(); self.angle_spin.setRange(-720, 720); self.angle_spin.setDecimals(3)
        self.angle_spin.setSuffix(" °"); self.angle_spin.setValue(180.0)
        trow.addWidget(self.target_check); trow.addWidget(self.angle_spin); trow.addStretch()
        layout.addLayout(trow)
        arow = QHBoxLayout()
        arow.addWidget(QLabel("   axis:"))
        self.axis_spins = []
        for name, v in zip("xyz", (1.0, 0.0, 0.0)):
            sp = QDoubleSpinBox(); sp.setRange(-1, 1); sp.setDecimals(3); sp.setSingleStep(0.1); sp.setValue(v)
            sp.setPrefix(name + " "); arow.addWidget(sp); self.axis_spins.append(sp)
        arow.addStretch()
        layout.addLayout(arow)

        self.var_list = QListWidget()
        self.var_list.setMaximumHeight(150)
        self.var_widgets = []
        if show_variables:
            head = QHBoxLayout()
            head.addWidget(QLabel("free variables:"))
            if variables_hint:
                h = QLabel(variables_hint); h.setStyleSheet(f"color: {MUTED_FG};"); h.setWordWrap(True)
                head.addWidget(h, stretch=1)
            b_all = QPushButton("all"); b_none = QPushButton("none")
            b_all.clicked.connect(lambda: self._check_all(True)); b_none.clicked.connect(lambda: self._check_all(False))
            head.addWidget(b_all); head.addWidget(b_none)
            layout.addLayout(head)
            layout.addWidget(self.var_list)
            self.var_widgets = [b_all, b_none, self.var_list]
        else:
            self.var_list.hide()
            hint = QLabel(variables_hint or "free variables: the parameters ticked 'free'")
            hint.setStyleSheet(f"color: {MUTED_FG};"); hint.setWordWrap(True)
            layout.addWidget(hint)

        brow = QHBoxLayout()
        self.button = QPushButton("Solve"); self.button.clicked.connect(self.toggle)
        brow.addWidget(self.button)
        self.settings_button = QPushButton("Settings…")
        self.settings_button.setToolTip("solver settings: evaluation limit, tolerances, goal weights")
        self.settings_button.clicked.connect(self.settingsRequested.emit)
        brow.addWidget(self.settings_button); brow.addStretch()
        layout.addLayout(brow)
        self.status = QLabel(""); self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.status)
        lock_heights(self)

    # ------------------------------------------------------------ variables ---
    def set_variables(self, items):
        """items: list of (key, label). Check states are kept for keys that stay; new keys are checked."""
        old = {}
        for i in range(self.var_list.count()):
            it = self.var_list.item(i)
            old[it.data(Qt.ItemDataRole.UserRole)] = it.checkState()
        self.var_list.blockSignals(True)
        self.var_list.clear()
        for key, label in items:
            it = QListWidgetItem(label)
            it.setData(Qt.ItemDataRole.UserRole, key)
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(old.get(key, Qt.CheckState.Checked))
            self.var_list.addItem(it)
        self.var_list.blockSignals(False)

    def free_keys(self):
        out = []
        for i in range(self.var_list.count()):
            it = self.var_list.item(i)
            if it.checkState() == Qt.CheckState.Checked:
                out.append(it.data(Qt.ItemDataRole.UserRole))
        return out

    def _check_all(self, on):
        for i in range(self.var_list.count()):
            self.var_list.item(i).setCheckState(Qt.CheckState.Checked if on else Qt.CheckState.Unchecked)

    def goals(self) -> Goals:
        orders = tuple(k for k, c in self.order_checks.items() if c.isChecked())
        target = None
        if self.target_check.isChecked():
            target = (tuple(sp.value() for sp in self.axis_spins), self.angle_spin.value())
        return Goals(orders, target, self.b1_check.isChecked())

    # --------------------------------------------------------------- process ---
    @property
    def solving(self):
        return self._process is not None

    def toggle(self):
        if self.solving:
            self.stop()
        else:
            self.start()

    def _set_inputs_enabled(self, on):
        widgets = list(self.order_checks.values()) + [self.b1_check, self.target_check, self.angle_spin,
                                                       self.settings_button] + self.axis_spins + self.var_widgets
        for w in widgets:
            w.setEnabled(on)

    def start(self):
        try:
            goals = self.goals()
            if goals.target is not None and sum(v * v for v in goals.target[0]) == 0:
                raise ValueError("the target rotation axis must not be zero")
            job = self.job_provider(goals, self.free_keys())
        except Exception as e:
            self._show(str(e), BAD)
            return
        job["goals"] = goals.to_dict()
        settings = self.settings_provider() if self.settings_provider else None
        job["settings"] = settings.to_dict() if settings is not None else None
        from psa.design.optimize import SolverSettings
        self._limit = (settings or SolverSettings()).max_nfev
        self._started = time.time()
        ctx = mp.get_context("spawn")
        self._queue, self._stop_event = ctx.Queue(), ctx.Event()
        self._process = ctx.Process(target=run_job, args=(job, self._queue, self._stop_event), daemon=True)
        self._process.start()
        self._stop_requested_at = None
        self._poll = QTimer(self); self._poll.setInterval(100); self._poll.timeout.connect(self._poll_queue)
        self._poll.start()
        self.button.setText("Stop")
        self._set_inputs_enabled(False)
        self._show("starting solver …")
        self.running.emit(True)

    def stop(self, wait=False):
        if not self.solving:
            return
        self._stop_event.set()
        if self._stop_requested_at is None:
            self._stop_requested_at = time.time()
        self._show("stopping …")
        if wait:
            self._process.join(5)
            if self._process.is_alive():
                self._process.terminate()
            self._finish()

    def _poll_queue(self):
        while True:
            try:
                msg = self._queue.get_nowait()
            except queue_module.Empty:
                break
            except (EOFError, OSError):
                msg = ("error", "the solver process ended unexpectedly")
            if msg[0] == "progress":
                self._show(f"solving … {msg[1]} / {self._limit} evaluations, {time.time() - self._started:.0f} s, "
                           f"best residual {msg[2]:.2e}  (click Stop to interrupt)")
            elif msg[0] == "done":
                self._finish()
                self._report(msg[2])
                self.solved.emit(msg[1], msg[2])
                return
            elif msg[0] == "error":
                self._finish()
                self._show(msg[1], BAD)
                return
        if self._stop_requested_at is not None and time.time() - self._stop_requested_at > 5:
            self._process.terminate()
            self._finish()
            self._show("solver terminated (no result)", BAD)
        elif self._process is not None and not self._process.is_alive() and self._queue.empty():
            self._finish()
            self._show("the solver process ended without a result", BAD)

    def _finish(self):
        if getattr(self, "_poll", None) is not None:
            self._poll.stop()
            self._poll = None
        if self._process is not None:
            self._process.join(1)
        self._process = None
        self.button.setText("Solve")
        self._set_inputs_enabled(True)
        self.running.emit(False)

    def _report(self, info):
        reason = info.get("stop_reason")
        if info["reached"]:
            state = "solved"
        else:
            state = {"user": "stopped – best so far", "limit": "evaluation limit reached – best so far",
                     "time": "time limit reached – best so far"}.get(reason, "not fully solved (converged)")
        values = ", ".join(f"{k} = {v:.2e}" for k, v in info["values"].items())
        start = ", ".join(f"{k} = {v:.2e}" for k, v in info.get("start_summary", {}).items())
        self._show(f"{state} ({info['nfev']} evaluations, {info.get('elapsed_s', 0):.1f} s): {values}\nbefore: {start}",
                   GOOD if info["reached"] else BAD)

    def _show(self, text, color=None):
        self.status.setStyleSheet(f"color: {color};" if color else "")
        self.status.setText(text)
