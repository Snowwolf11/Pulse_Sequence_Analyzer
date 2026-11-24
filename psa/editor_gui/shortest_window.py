"""Small window of the GRAPE tab: the shortest robust pulse (psa.design.shortest). The goal comes from the tab
(gate goal, offset orders, B1 1st order; on resonance, the editor's dt and rf amplitude); the window holds the
search settings, runs the search in a background process and shows its log. Every accepted (shorter) pulse
goes to the tab's plots (pulseFound); the shortest one stays there as the tab's result."""
from __future__ import annotations

import multiprocessing as mp
import os
import queue as queue_module

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (QCheckBox, QDialog, QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel,
                               QPlainTextEdit, QPushButton, QSpinBox, QVBoxLayout)

from psa.design.pulse import PulseSequence
from psa.design.shortest import ShortestSpec, search_process
from psa.editor_gui.widgets import BAD, GOOD, MUTED_FG

INTRO = ("Shortest pulse found that does the gate exactly and cancels the offset error through the ticked orders "
         "(and first-order B1), at most 100 % rf amplitude. GRAPE (local mode) solves the goals at a start duration, "
         "then shortens the pulse step by step (10 %, then bisection) from the best solution and from fresh starts "
         "(random, CORPSE, 2H / 4H / 6H); every accepted pulse is checked independently. A local search: the "
         "shortest pulse found, not proven optimal (in the reliability check, 3 seeds and start durations gave the "
         "same minimum in all 12 test cases).")


class ShortestWindow(QDialog):
    pulseFound = Signal(object)              # PulseSequence: every accepted pulse and the final one
    runningChanged = Signal(bool)

    def __init__(self, spec_provider, start_provider=None, parent=None):
        """spec_provider() -> (ShortestSpec, description) or raises ValueError; start_provider() -> the tab's
        start pulse (or None)."""
        super().__init__(parent)
        self.setWindowTitle("Shortest robust pulse")
        self.setWindowFlag(Qt.WindowType.Window, True)
        self.resize(620, 560)
        self.spec_provider, self.start_provider = spec_provider, start_provider
        self._proc = None
        self.result, self.verification, self.last_result = None, None, None
        lay = QVBoxLayout(self)
        intro = QLabel(INTRO); intro.setWordWrap(True); intro.setStyleSheet(f"color: {MUTED_FG};")
        lay.addWidget(intro)
        f = QFormLayout()
        self.goal_label = QLabel(""); self.goal_label.setWordWrap(True)
        self.goal_label.setToolTip("set in the GRAPE tab: gate goal, offset orders (1st, 1st + 2nd or 1st – 3rd) and "
                                   "B1 1st order; the points of the tab are not used (the orders expand around 0 kHz)")
        f.addRow("goal (from the tab):", self.goal_label)
        self.start_spin = QDoubleSpinBox(); self.start_spin.setRange(1, 1e5); self.start_spin.setDecimals(1)
        self.start_spin.setSuffix(" µs"); self.start_spin.setValue(300.0)
        self.start_spin.setToolTip("first duration; should be long enough for a solution (if none is found, the "
                                   "search lengthens it by 25 % up to 4 times). The result does not depend much on it")
        f.addRow("start duration:", self.start_spin)
        self.time_spin = QDoubleSpinBox(); self.time_spin.setRange(0.5, 1440); self.time_spin.setDecimals(1)
        self.time_spin.setSuffix(" min"); self.time_spin.setValue(30.0)
        self.time_spin.setToolTip("the search ends after this time with the shortest pulse so far (checked between "
                                  "durations). In the reliability check a search took 4–36 min on one core per "
                                  "candidate (the most for 3rd order + B1)")
        f.addRow("time limit:", self.time_spin)
        row = QHBoxLayout()
        self.workers_spin = QSpinBox(); self.workers_spin.setRange(1, 256)
        self.workers_spin.setValue(max(1, min(10, (os.cpu_count() or 2) - 1)))
        self.workers_spin.setToolTip("processes for the candidates of one duration (up to 10 at a time)")
        self.seed_spin = QSpinBox(); self.seed_spin.setRange(0, 10 ** 6); self.seed_spin.setPrefix("seed ")
        self.seed_spin.setToolTip("changes every random start and perturbation: repeat with another seed (and start "
                                  "duration) to check the result")
        row.addWidget(self.workers_spin); row.addWidget(self.seed_spin); row.addStretch()
        f.addRow("processes:", row)
        self.use_start_check = QCheckBox("also try the tab's start pulse (squeezed to every duration)")
        f.addRow("", self.use_start_check)
        lay.addLayout(f)
        brow = QHBoxLayout()
        self.run_button = QPushButton("Search"); self.run_button.clicked.connect(self.toggle)
        self.save_button = QPushButton("Save…"); self.save_button.clicked.connect(self.save)
        self.save_button.setEnabled(False)
        brow.addWidget(self.run_button); brow.addWidget(self.save_button); brow.addStretch()
        lay.addLayout(brow)
        self.status = QLabel(""); self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.status)
        self.log = QPlainTextEdit(); self.log.setReadOnly(True)
        self.log.setToolTip("one line per duration tried: ok / fail, solved candidates, the best start, time")
        lay.addWidget(self.log, 1)
        self.input_widgets = [self.start_spin, self.time_spin, self.workers_spin, self.seed_spin, self.use_start_check]
        self._timer = QTimer(self); self._timer.setInterval(200); self._timer.timeout.connect(self.poll)
        self.refresh_goal()

    # ------------------------------------------------------------- inputs ---
    def refresh_goal(self):
        """Show the goal of the tab (or why it cannot be used); returns the spec or None."""
        try:
            spec, text = self.spec_provider()
        except ValueError as e:
            self.goal_label.setText(str(e)); self.goal_label.setStyleSheet(f"color: {BAD};")
            return None
        self.goal_label.setText(f"{text}; lower bound {spec.lower_bound_s() * 1e6:.1f} µs")
        self.goal_label.setStyleSheet("")
        return spec

    def showEvent(self, event):
        super().showEvent(event)
        if not self.running:
            self.refresh_goal()

    # ---------------------------------------------------------------- run ---
    @property
    def running(self):
        return self._proc is not None

    def toggle(self):
        self.stop() if self.running else self.start()

    def start(self):
        spec = self.refresh_goal()
        if spec is None:
            self._show("the goal of the GRAPE tab cannot be used (see above)", BAD)
            return
        start_arr = None
        if self.use_start_check.isChecked():
            ps = self.start_provider() if self.start_provider else None
            if ps is None:
                self._show("the GRAPE tab has no start pulse", BAD)
                return
            start_arr = np.column_stack([ps.amplitude, ps.phase]).astype(float)
        self.spec = spec
        self.result, self.verification, self.last_result = None, None, None
        self.save_button.setEnabled(False)
        self.log.clear()
        ctx = mp.get_context("spawn")
        self._queue, self._stop = ctx.Queue(), ctx.Event()
        kwargs = dict(time_limit_s=self.time_spin.value() * 60.0, seed=self.seed_spin.value())
        self._proc = ctx.Process(target=search_process, daemon=False,          # it starts its own pool
                                 args=(dict(spec.__dict__), self.start_spin.value(), kwargs, self._queue, self._stop,
                                       self.workers_spin.value(), start_arr))
        self._proc.start()
        self._stop_asked = False
        for w in self.input_widgets:
            w.setEnabled(False)
        self.run_button.setText("Stop")
        self._show("searching … (the GRAPE tab shows every accepted pulse)")
        self._timer.start()
        self.runningChanged.emit(True)

    def stop(self, wait=False):
        if not self.running:
            return
        self._stop.set()
        self._stop_asked = True
        self._show("stopping … (keeps the shortest pulse so far)")
        if wait:
            self._proc.join(15)
            self.poll()
            if self.running:
                self._proc.terminate()
                self._end("stopped")

    def poll(self):
        if not self.running:
            return
        while True:
            try:
                msg = self._queue.get_nowait()
            except queue_module.Empty:
                break
            except (EOFError, OSError):
                self._end("the search process ended unexpectedly", BAD)
                return
            kind = msg[0]
            if kind == "log":
                self.log.appendPlainText(msg[1])
            elif kind == "best":
                self._accept(msg[1], msg[2], msg[3])
                v = msg[3]
                self._show(f"searching … shortest so far {v['duration [µs]']:.1f} µs (N = {v['N']})")
            elif kind == "done":
                res, amp, phase = msg[1], msg[2], msg[3]
                self.last_result = res
                if amp is not None:
                    self._accept(amp, phase, res["verification"])
                self._end(self._summary(res), GOOD if amp is not None else BAD)
                return
            elif kind == "error":
                self.log.appendPlainText(msg[1])
                self._end("the search failed: " + msg[1].splitlines()[0], BAD)
                return
        if self.running and not self._proc.is_alive() and self._queue.empty():
            self._end("the search process ended without a result", BAD)

    def _accept(self, amp, phase, v):
        ps = PulseSequence(np.asarray(amp, float), np.asarray(phase, float), self.spec.dt, self.spec.umax)
        self.result, self.verification = ps, v
        self.save_button.setEnabled(True)
        self.pulseFound.emit(ps)

    def _summary(self, res):
        v = res["verification"]
        head = {"stopped": "stopped", "time limit": "time limit reached"}.get(res["message"], "finished")
        if v is None:
            return f"{head}: no solution found ({res['message']}, {res['elapsed_s'] / 60:.1f} min)"
        orders = ", ".join(f"{k} {v[k]:.1e}" for k in ("A1", "A2", "A3")[:self.spec.order])
        return (f"{head}: shortest pulse found {v['duration [µs]']:.1f} µs (N = {v['N']}); lower bound "
                f"{res['lower bound [µs]']:.1f} µs; {res['message']}; {res['elapsed_s'] / 60:.1f} min\n"
                f"check: gate error {v['gate error']:.1e}" + (f", {orders}" if orders else "")
                + (f", B1 {v['B1']:.1e}" if self.spec.b1 else "") + f", max. amplitude {v['max amplitude [%]']:.1f} %"
                + (" – passed" if v["ok"] else " – NOT passed"))

    def _end(self, text, color=None):
        self._timer.stop()
        if self._proc is not None:
            self._proc.join(1)
            if self._proc.is_alive():
                self._proc.terminate()
        self._proc = None
        for w in self.input_widgets:
            w.setEnabled(True)
        self.run_button.setText("Search")
        self._show(text, color)
        self.runningChanged.emit(False)

    def save(self):
        if self.result is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save pulse", f"shortest_{self.result.duration * 1e6:.0f}us.bruker",
                                              "Bruker shape (*.bruker);;CSV (*.csv)")
        if path:
            self.result.save(path)
            self._show(self.status.text() + f"\nsaved: {path}", self.status.property("color_") or None)

    def _show(self, text, color=None):
        self.status.setStyleSheet(f"color: {color};" if color else "")
        self.status.setProperty("color_", color)
        self.status.setText(text)

    def closeEvent(self, event):              # closing the window hides it; a running search goes on
        event.ignore()
        self.hide()
