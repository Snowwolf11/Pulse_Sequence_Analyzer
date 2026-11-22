"""GUI tests of the pulse-sequence editor (PySide6 + pyqtgraph, run offscreen).

Run with ``python -m pytest psa/tests`` or ``python psa/tests/test_editor_gui.py``.
Real mouse interaction is simulated with QTest on the plot widgets.
"""
import os
import sys
import tempfile

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from PySide6.QtCore import QPoint, QPointF, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from psa.design.pulse import PulseSequence
from psa.design.conversion import ps_to_curve
from psa.design import metrics

APP = QApplication.instance() or QApplication(sys.argv)


def process(ms=0):
    if ms:
        QTest.qWait(ms)
    APP.processEvents()


def make_window(main=None):
    from psa.editor_gui.window import PulseEditorWindow
    w = PulseEditorWindow(main_window=main, dt=0.5e-6, umax=1e4, save_path=os.path.join(tempfile.mkdtemp(), "ed.csv"))
    w.resize(1500, 900)
    w.show()
    process(50)
    return w


def wait_panel(w):
    w.panel._timer.stop()
    w.panel.recompute()


def view_to_widget(plot, x, y):
    vb = plot.getPlotItem().vb
    scene_pt = vb.mapViewToScene(QPointF(x, y))
    p = plot.mapFromScene(scene_pt)
    return QPoint(int(p.x()), int(p.y()))


def drag(plot, start, end, steps=8):
    vp = plot.viewport()
    QTest.mousePress(vp, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, start)
    for i in range(1, steps + 1):
        QTest.mouseMove(vp, QPoint(start.x() + (end.x() - start.x()) * i // steps, start.y() + (end.y() - start.y()) * i // steps))
        process(5)
    QTest.mouseRelease(vp, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, end)
    process(20)


def double_click(plot, pos):
    """The event sequence of a real double-click: press, release, double-click, release."""
    vp, L, NM = plot.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier
    QTest.mousePress(vp, L, NM, pos, 500)
    QTest.mouseRelease(vp, L, NM, pos, 10)
    QTest.mouseDClick(vp, L, NM, pos, 10)
    QTest.mouseRelease(vp, L, NM, pos, 10)
    process(20)


def test_text_tab_and_panel():
    w = make_window()
    w.text_tab.set_text("R(100, 0, 25); R(100, 90, 50);")
    wait_panel(w)
    ps = w.panel.pulse
    assert ps.N == 150
    g = metrics.gate(ps)
    assert w.panel.labels["gate"].text().startswith(f"{g['angle_deg']:.4f}")
    w.text_tab.editor.setPlainText("R(100, 0, 0.3);")                     # not a multiple of dt
    w.text_tab.parse()
    assert w.text_tab.error.text() != ""
    w.close()


def test_waveform_tab_keyframes_and_duration():
    w = make_window()
    w.show_tab(w.waveform_tab)
    tab = w.waveform_tab
    tab.duration_spin.setValue(0.05); tab._duration_edited()                # 50 us -> 100 sub-pulses
    assert tab.model.N == 100
    tab.add_keyframe("amplitude", 50, 40.0)
    tab.add_keyframe("phase", 99, 720.0)
    ps = tab.current_pulse()
    assert ps.N == 100 and abs(ps.amplitude[50] - 40.0) < 1e-12 and ps.phase[-1] == 720.0
    assert abs(ps.amplitude[25] - 70.0) < 1e-12                             # linear 100 -> 40
    tab._select("amplitude", tab.model.amplitude.find(0))
    tab.sel_mode.setCurrentIndex(1)                                        # hold until next keyframe
    assert np.allclose(tab.current_pulse().amplitude[:50], 100.0)
    wait_panel(w)
    assert w.panel.pulse.N == 100
    w.close()


def test_mouse_drag_double_click_and_delete():
    w = make_window()
    w.show_tab(w.waveform_tab)
    tab = w.waveform_tab
    plot = tab.plots["amplitude"]
    tab.add_keyframe("amplitude", 100, 100.0)                              # keyframe at 50 us
    plot.getPlotItem().vb.setRange(xRange=(0, 100), yRange=(0, 105), padding=0)
    process(20)
    start = view_to_widget(plot, 50.0, 100.0)
    end = view_to_widget(plot, 70.3, 30.0)
    drag(plot, start, end)
    k = tab.model.amplitude.keyframes[tab.model.amplitude.find(141) if tab.model.amplitude.find(141) is not None else -1]
    assert k.index == 141, [(kk.index, kk.value) for kk in tab.model.amplitude.keyframes]   # 70.3 us -> sub-pulse 141
    assert abs(k.value - 30.0) < 2.0, k.value
    # double-click adds a keyframe
    n_before = len(tab.model.amplitude.keyframes)
    double_click(plot, view_to_widget(plot, 20.0, 60.0))
    process(20)
    kfs = [(kk.index, round(kk.value, 1)) for kk in tab.model.amplitude.keyframes]
    assert len(tab.model.amplitude.keyframes) == n_before + 1, kfs
    assert tab.model.amplitude.find(40) is not None, kfs
    # right-click deletes it again
    QTest.mouseClick(plot.viewport(), Qt.MouseButton.RightButton, Qt.KeyboardModifier.NoModifier, view_to_widget(plot, 20.0, tab.model.amplitude.keyframes[tab.model.amplitude.find(40)].value))
    process(20)
    assert tab.model.amplitude.find(40) is None, [(kk.index, round(kk.value, 1)) for kk in tab.model.amplitude.keyframes]
    w.close()


def test_waveform_axis_ranges():
    w = make_window()
    w.show_tab(w.waveform_tab)
    tab = w.waveform_tab
    tab.duration_spin.setValue(0.2); tab._duration_edited()
    tab.set_ranges(time=(20, 80), amplitude=(10, 60), phase=(-90, 270))
    (x0, x1), (a0, a1) = tab.plots["amplitude"].getPlotItem().vb.viewRange()
    (px0, px1), (p0, p1) = tab.plots["phase"].getPlotItem().vb.viewRange()
    assert np.allclose([x0, x1], [20, 80], atol=1e-6) and np.allclose([px0, px1], [20, 80], atol=1e-6)   # time shared
    assert np.allclose([a0, a1], [10, 60], atol=1e-6) and np.allclose([p0, p1], [-90, 270], atol=1e-6)
    # zooming the view (as with the trackpad) updates the fields
    tab.plots["phase"].getPlotItem().vb.setYRange(0, 720, padding=0); process(20)
    assert abs(tab.range_spins["phase"][1].value() - 720) < 1e-6
    tab.auto_range()
    (x0, x1), _ = tab.plots["amplitude"].getPlotItem().vb.viewRange()
    assert x0 <= 0 and x1 >= 200
    w.close()


def test_waveform_import_from_text_and_metadata_change():
    w = make_window()
    w.text_tab.set_text("R(80, 0, 20); PR(100, 0, 270, 30); R(50, 180, 10);")
    w.show_tab(w.waveform_tab)
    w.waveform_tab.amp_tol.setValue(0); w.waveform_tab.phase_tol.setValue(0)
    w._waveform_from_text()
    a, b = w.text_tab.current_pulse(), w.waveform_tab.current_pulse()
    assert np.allclose(a.amplitude, b.amplitude) and np.allclose((a.phase - b.phase + 180) % 360 - 180, 0, atol=1e-9)
    # change dt in the panel: waveform keeps its duration, text tab re-parses with the new dt
    w.panel.dt_spin.setValue(1.0); w.panel._metadata_edited()
    assert w.waveform_tab.model.N == 60 and abs(w.waveform_tab.model.duration - 60e-6) < 1e-12
    assert w.text_tab.current_pulse().N == 60
    w.close()


def test_families_tab():
    from psa.design.families import family_by_name
    from psa.design.conversion import curve_to_ps
    w = make_window()
    tab = w.families_tab
    w.show_tab(tab)
    for name in ("Helix", "Torus knot", "Custom (formula)", "2H (closed two-helix loop)", "4H (0° UR, 2nd-order robust)"):
        tab.family_box.setCurrentText(name)
        tab._timer.stop()
        ps = tab.compute()
        assert ps is not None and ps.N > 0, (name, tab.error.text())
        wait_panel(w)
        assert w.panel.pulse is ps
    # parameter change through the slider control changes the pulse
    tab.family_box.setCurrentText("Helix"); tab._timer.stop(); n1 = tab.compute().N
    tab.controls["turns"].set_value(4.0); tab.controls["turns"].changed.emit(); tab._timer.stop()
    assert tab.compute().N > n1
    # same result as the backend
    ref, _ = curve_to_ps(family_by_name("Helix").points(tab.values()), 0.5e-6, 1e4)
    assert np.allclose(tab.pulse.amplitude, ref.amplitude)
    # fixed duration
    tab.dur_fixed.setChecked(True); tab.dur_spin.setValue(1.0); tab._timer.stop()
    assert abs(tab.compute().duration - 1e-3) < 1e-12
    # a bad formula shows an error and keeps the last pulse
    tab.family_box.setCurrentText("Custom (formula)"); tab._timer.stop(); tab.compute()
    tab.expr_edits[0].setText("cos(u"); tab._timer.stop()
    assert tab.compute() is None and "syntax" in tab.error.text()
    # design round trip and hand-over to the waveform tab
    tab.expr_edits[0].setText("cos(a*u)"); tab.family_box.setCurrentText("4H (0° UR, 2nd-order robust)")
    tab.controls["tau"].set_value(70.0); tab._timer.stop(); ps = tab.compute()
    d = tab.to_dict()
    tab.family_box.setCurrentText("Helix"); tab.load_design(d)
    assert np.allclose(tab.pulse.phase, ps.phase)
    w.waveform_tab.amp_tol.setValue(0); w.waveform_tab.phase_tol.setValue(0)
    w._edit_in_waveform(ps)
    assert w.current_tab() is w.waveform_tab
    assert np.allclose(w.waveform_tab.current_pulse().amplitude, ps.amplitude)
    w.close()


def test_layout_is_scrollable_and_fits_screen():
    from PySide6.QtWidgets import QScrollArea
    from PySide6.QtGui import QGuiApplication
    w = make_window()
    w._fit_to_screen()
    avail = QGuiApplication.primaryScreen().availableGeometry()
    # the offscreen test platform reports an 800 px wide screen, narrower than the window's minimum width
    assert w.height() <= avail.height() and w.width() <= max(avail.width(), w.minimumSizeHint().width())
    assert w.minimumSizeHint().width() <= 1280                             # fits a 13" laptop screen
    for i in range(w.tabs.count()):
        assert isinstance(w.tabs.widget(i), QScrollArea) and w.tabs.widget(i).widget() is w.all_tabs()[i]
    assert isinstance(w.panel.parentWidget().parentWidget(), QScrollArea)
    # a small window shows a scroll bar instead of squeezing the panel
    w.resize(900, 450); process(50)
    area = w.panel.parentWidget().parentWidget()
    assert area.verticalScrollBar().maximum() > 0
    w.close()


def test_save_format_selector():
    w = make_window()
    d = tempfile.mkdtemp()
    w.panel.path_edit.setText(os.path.join(d, "pulse.csv")); w.panel._sync_format_from_path()
    assert w.panel.format_box.currentIndex() == 0
    w.panel.format_box.setCurrentIndex(1)                                  # Bruker
    assert w.panel.path_edit.text().endswith("pulse.bruker")
    w.text_tab.set_text("R(100, 0, 25);"); wait_panel(w)
    path = w.panel.save()
    assert path.endswith(".bruker") and open(path).read().startswith("##TITLE")
    w.panel.path_edit.setText(os.path.join(d, "noext"))                    # no extension: selected format is used
    assert w.panel.save().endswith("noext.bruker")
    w.close()


def test_families_live_switch():
    w = make_window()
    tab = w.families_tab
    w.show_tab(tab)
    tab.family_box.setCurrentText("Lissajous knot"); tab._timer.stop(); tab._update()
    assert tab.pulse is not None
    tab.live_check.setChecked(False)
    assert tab.pulse is None and "not calculated" in tab.info.text()
    wait_panel(w)
    assert w.panel.pulse is None and w.panel.labels["info"].text() == "not calculated"
    import time
    tab.controls["nz"].set_value(5); tab.controls["nz"].changed.emit(); tab._timer.stop()
    t0 = time.perf_counter(); tab._update(); dt_curve = time.perf_counter() - t0
    assert tab.pulse is None                                               # only the curve was drawn
    assert dt_curve < 1.0, dt_curve
    tab.calc_button.click()
    assert tab.pulse is not None
    wait_panel(w)
    assert w.panel.pulse is tab.pulse
    tab.live_check.setChecked(True)
    w.close()


def test_builder_tab():
    from psa.design.builder import Segment
    w = make_window()
    tab = w.builder_tab
    w.show_tab(tab)
    tab._timer.stop(); ps = tab.recompute()
    assert ps is not None and ps.N > 0
    wait_panel(w)
    assert w.panel.pulse is ps
    # add segments through the buttons and edit parameters through the editor widgets
    n0 = len(tab.design.segments)
    tab.add_segment("helix"); tab.add_segment("line"); tab.add_segment("bezier")
    assert len(tab.design.segments) == n0 + 3 and tab.seg_list.count() == n0 + 3
    tab.seg_list.setCurrentRow(0)
    spin, free = tab.editors["angle"]
    spin.setValue(250.0)
    assert tab.design.segments[0].params["angle"] == 250.0
    tab.remove_segment(); assert len(tab.design.segments) == n0 + 2
    tab.move_segment(1)
    # solver through the GUI: three arcs, twists and angles free -> closed curve
    rmin = tab.design.r_min_us
    tab.design.segments = [Segment("arc", {"radius": 1.3 * rmin, "angle": 90}),
                           Segment("arc", {"radius": 1.3 * rmin, "angle": 200}, twist=150, free=["angle", "twist"]),
                           Segment("arc", {"radius": 1.3 * rmin, "angle": 120}, twist=150, free=["angle", "twist"])]
    tab._refresh_list(select=0)
    tab.solve_button.click()                                                # runs in the background
    assert tab.solving and tab.solve_button.text() == "Stop" and tab.solve_button.isEnabled()
    assert not tab.seg_list.isEnabled()                                     # editing locked while solving
    wait_until(lambda: not tab.solving, 60)
    assert tab.solve_label.text().startswith("solved"), tab.solve_label.text()
    assert tab.seg_list.isEnabled() and tab.solve_button.text() == "Solve"
    tab._timer.stop(); tab.recompute()
    assert tab.design.robustness()["A1_norm"] < 1e-8                         # closed (below the "solved" threshold)
    assert np.linalg.norm(tab.points[-1]) < 1e-4 * tab.design.duration()    # drawn curve (0.25 µs chords)
    # design round trip
    d = tab.design.to_dict()
    tab.load_design(d)
    assert tab.design.robustness()["A1_norm"] < 1e-8
    w.close()


def wait_until(cond, timeout_s):
    import time
    t0 = time.time()
    while not cond():
        process(20)
        assert time.time() - t0 < timeout_s, "timeout"


def test_builder_solver_stop_keeps_gui_responsive():
    """The solver runs in a separate process: the GUI keeps running and Stop ends it with a result."""
    import time
    from psa.design.builder import Segment
    w = make_window()
    tab = w.builder_tab
    w.show_tab(tab)
    rmin = tab.design.r_min_us
    tab.design.segments = [Segment("arc", {"radius": 1.3 * rmin, "angle": 100 + 37 * k}, twist=60 * k,
                                   free=(["angle", "twist", "radius"] if k else [])) for k in range(7)]
    tab._refresh_list(select=0)
    for k in (1, 2, 3):
        tab.order_checks[k].setChecked(True)
    t0, ticks = time.time(), 0
    while time.time() - t0 < 1.0:                                          # idle event-loop rate for comparison
        process(10); ticks += 1
    idle = ticks / (time.time() - t0)
    tab.solve_button.click()
    assert tab.solving
    t0, ticks = time.time(), 0
    while time.time() - t0 < 1.5 and tab.solving:                          # GUI event loop keeps running
        process(10); ticks += 1
    rate = ticks / (time.time() - t0)
    assert rate > 0.5 * idle, (rate, idle)
    if tab.solving:
        tab.solve_button.click()                                            # Stop
        wait_until(lambda: not tab.solving, 30)
    text = tab.solve_label.text()
    assert text.startswith(("stopped", "solved", "not fully")), text
    assert tab.seg_list.isEnabled() and tab.solve_button.text() == "Solve"
    w.close()


def test_waveform_solver():
    """Waveform tab: the solver adjusts keyframe values for a target gate; plots keep their size."""
    from psa.design import metrics
    w = make_window()
    tab = w.waveform_tab
    w.show_tab(tab)
    tab.duration_spin.setValue(0.12); tab._duration_edited()            # 240 sub-pulses
    for name, vals in (("amplitude", (80.0, 60.0, 50.0)), ("phase", (0.0, 60.0, 30.0))):
        for idx, v in zip((0, 120, 239), vals):
            tab.add_keyframe(name, idx, v, "hold")
    keys = tab.solver_box.free_keys()
    assert ("amplitude", 0) in keys and ("phase", 120) in keys
    # do not let the solver change the last phase keyframe (one sub-pulse)
    for i in range(tab.solver_box.var_list.count()):
        it = tab.solver_box.var_list.item(i)
        if it.data(Qt.ItemDataRole.UserRole) == ("phase", 239):
            it.setCheckState(Qt.CheckState.Unchecked)
    box = tab.solver_box
    box.order_checks[1].setChecked(False)
    box.target_check.setChecked(True); box.angle_spin.setValue(90.0)
    for sp, v in zip(box.axis_spins, (0, 1, 0)):
        sp.setValue(v)
    box.button.click()
    assert box.solving and not tab.duration_spin.isEnabled()
    wait_until(lambda: not box.solving, 120)
    assert box.status.text().startswith("solved"), box.status.text()
    g = metrics.gate(tab.current_pulse())
    assert abs(g["angle_deg"] - 90) < 1e-3 and np.allclose(g["axis"], [0, 1, 0], atol=1e-4)
    assert tab.model.phase.keyframes[tab.model.phase.find(239)].value == 30.0
    assert tab.duration_spin.isEnabled()
    # the plots fill the visible area; the solver box lies below
    area = w.tabs.currentWidget()
    w.resize(1300, 800); process(50)
    assert tab.plots["phase"].geometry().bottom() <= area.viewport().height() + 5
    assert tab.solver_box.geometry().top() >= tab.plots["phase"].geometry().bottom()
    w.close()


def test_families_solver():
    w = make_window()
    tab = w.families_tab
    w.show_tab(tab)
    tab.family_box.setCurrentText("Custom (formula)")
    tab.expr_edits[1].setText("sin(b*u) + a*u")
    tab.controls["a"].set_value(0.3); tab.controls["b"].set_value(1.0)
    tab._timer.stop(); tab.compute()
    names = [tab.solver_box.var_list.item(i).data(Qt.ItemDataRole.UserRole) for i in range(tab.solver_box.var_list.count())]
    assert names == ["u0", "u1", "a", "b", "c", "d"]
    tab.solver_box._check_all(False)
    tab.solver_box.var_list.item(2).setCheckState(Qt.CheckState.Checked)     # only "a"
    duration = tab.pulse.duration
    tab.solver_box.button.click()
    wait_until(lambda: not tab.solver_box.solving, 120)
    assert abs(tab.controls["a"].value()) < 1e-3, tab.solver_box.status.text()
    assert abs(tab.pulse.duration - duration) < 1e-12                      # duration kept
    w.close()


def test_builder_bezier_drag():
    from psa.design.builder import Segment
    w = make_window()
    tab = w.builder_tab
    w.show_tab(tab)
    tab.design.segments = [Segment("arc", {"radius": 25, "angle": 60}),
                           Segment("bezier", {"p1": 10, "p2x": 8, "p2y": 2, "p2z": 20, "p3x": 20, "p3y": -4, "p3z": 30})]
    tab._refresh_list(select=1); tab._timer.stop(); tab.recompute(); process(30)
    handles = tab._control_points() * 1e6
    plot = tab.proj_plots["xy"]
    plot.getPlotItem().vb.autoRange(); process(30)
    before = dict(tab.design.segments[1].params)
    start = view_to_widget(plot, handles[2, 0], handles[2, 1])                 # P3 in the xy projection
    end = view_to_widget(plot, handles[2, 0] + 5.0, handles[2, 1] - 3.0)
    drag(plot, start, end)
    after = tab.design.segments[1].params
    moved = tab._control_points()[2] * 1e6
    assert abs(moved[0] - (handles[2, 0] + 5.0)) < 0.6 and abs(moved[1] - (handles[2, 1] - 3.0)) < 0.6, (moved, handles[2])
    assert abs(moved[2] - handles[2, 2]) < 1e-6                                # z unchanged in the xy projection
    assert after["p2x"] == before["p2x"] and after["p1"] == before["p1"]       # other points untouched
    w.close()


def test_import_tab():
    from psa.design.families import family_by_name
    from psa.design.conversion import curve_to_ps
    w = make_window()
    tab = w.import_tab
    w.show_tab(tab)
    pts = family_by_name("Torus knot").points({})
    path = os.path.join(tempfile.mkdtemp(), "knot.csv")
    np.savetxt(path, pts, delimiter=",", header="x,y,z")
    tab.path_edit.setText(path); tab.format_box.setCurrentIndex(0)
    ps = tab.load()
    ref, _ = curve_to_ps(pts, 0.5e-6, 1e4)
    assert ps is not None and np.allclose(ps.amplitude, ref.amplitude) and np.allclose(ps.phase, ref.phase)
    wait_panel(w)
    assert w.panel.pulse is ps
    tab.dur_fixed.setChecked(True); tab.dur_spin.setValue(2.0)
    assert abs(tab.convert().duration - 2e-3) < 1e-12
    tab.path_edit.setText(path + ".missing")
    assert tab.load() is None and "could not read" in tab.error.text()
    w.close()


def test_save_and_send_to_main_window():
    from psa.GUI_pyside6 import PulseSequenceAnalyzerQt
    main = PulseSequenceAnalyzerQt()
    main.show(); process(20)
    main.inpo_fact_spin.setValue(1)                        # same sampling as the editor preview
    w = make_window(main)
    w.text_tab.set_text("R(100, 0, 25); PR(100, 0, 360, 50); R(100, 180, 25);")
    wait_panel(w)
    for ext in (".csv", ".bruker"):
        path = os.path.join(tempfile.mkdtemp(), "pulse" + ext)
        w.panel.path_edit.setText(path)
        w.panel.save_and_send()
        process(50)
        assert os.path.exists(path)
        back = PulseSequence.load(path)
        assert np.allclose(back.amplitude, w.panel.pulse.amplitude)
        assert main.pulse_sequence_edit.text() == path
        assert main.latest_result is not None
        # the analyzer's curve equals the editor's curve
        assert np.allclose(main.latest_result.CM[-1], ps_to_curve(w.panel.pulse)[-1], atol=1e-12)
    w.close(); main.close()


def test_help_window():
    w = make_window()
    w.help_button.click(); process(20)
    from PySide6.QtWidgets import QTextEdit
    assert w._help.isVisible() and "Curve builder tab" in w._help.findChild(QTextEdit).toPlainText()
    from psa.documentation import DOCUMENTATION_STR
    assert "Pulse Sequence Editor" in DOCUMENTATION_STR and "Solver" in DOCUMENTATION_STR
    w.close()


def test_solver_settings_shared_and_used():
    """Settings → Solver of the analyzer edit the settings the editor's solver uses."""
    from psa.GUI_pyside6 import PulseSequenceAnalyzerQt, SettingsDialog
    from psa.design.builder import Segment
    main = PulseSequenceAnalyzerQt()
    main.open_editor(); process(50)
    ed = main._editor_window
    assert ed.solver_settings is main.solver_settings
    dlg = SettingsDialog(main, initial_tab="Solver")
    assert dlg.tabs.count() == 2 and dlg.tabs.currentWidget() is dlg.solver_widget
    dlg.solver_widget.max_spin.setValue(10)
    dlg.solver_widget.time_spin.setValue(120.0)
    dlg.solver_widget.solved_check.setChecked(True)
    dlg.solver_widget.weights["weight_order1"].setValue(3.0)
    dlg.solver_widget.sci["solved_tolerance"].setText("1e-6")
    dlg.save()
    assert main.solver_settings.max_nfev == 10 and main.solver_settings.weight_order1 == 3.0
    assert main.solver_settings.time_limit_s == 120.0 and main.solver_settings.stop_when_solved
    assert main.solver_settings.solved_tolerance == 1e-6
    # an invalid entry is refused and keeps the dialog on the solver tab
    dlg2 = SettingsDialog(main)
    dlg2.solver_widget.sci["xtol"].setText("abc")
    dlg2.save()
    assert dlg2.solver_widget.error.text() and main.solver_settings.xtol == 1e-12
    # the job sent to the solver process carries the settings: 10 evaluations -> stops early
    tab = ed.builder_tab
    ed.show_tab(tab)
    rmin = tab.design.r_min_us
    tab.design.segments = [Segment("arc", {"radius": 1.3 * rmin, "angle": 100 + 37 * k}, twist=60 * k,
                                   free=(["angle", "twist", "radius"] if k else [])) for k in range(7)]
    tab._refresh_list(select=0)
    for k in (1, 2, 3):
        tab.order_checks[k].setChecked(True)
    tab.solve_button.click()
    wait_until(lambda: not tab.solving, 60)
    import re
    n = int(re.search(r"\((\d+) evaluations", tab.solve_label.text()).group(1))
    assert n == 10 and tab.solve_label.text().startswith("evaluation limit reached"), tab.solve_label.text()
    # Reset Analyzer resets the shared object in place
    main.reset_to_defaults()
    assert ed.solver_settings is main.solver_settings and main.solver_settings.max_nfev == 2000
    ed.close(); main.close()


def test_solver_settings_button_opens_dialog():
    from PySide6.QtWidgets import QDialog
    w = make_window()                                     # no analyzer: stand-alone dialog
    opened = []
    orig = QDialog.exec
    QDialog.exec = lambda self: opened.append(type(self).__name__) or 0
    try:
        w.waveform_tab.solver_box.settings_button.click()
    finally:
        QDialog.exec = orig
    assert opened == ["SolverSettingsDialog"]
    w.close()


def test_open_editor_from_main_window():
    from psa.GUI_pyside6 import PulseSequenceAnalyzerQt
    main = PulseSequenceAnalyzerQt()
    main.time_per_pulse_spin.setValue(1.0); main.pulse_amplitude_spin.setValue(20)
    main.open_editor(); process(50)
    ed = main._editor_window
    assert ed.isVisible() and abs(ed.panel.dt - 1e-6) < 1e-15 and abs(ed.panel.umax - 2e4) < 1e-9
    ed.close(); main.close()


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t(); print(f"[PASS] {t.__name__}")
        except Exception as e:
            failed += 1; print(f"[FAIL] {t.__name__}: {type(e).__name__}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
