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
    for name in ("Helix", "Torus knot", "Custom (formula)", "2H (closed two-helix loop)", "4H (0° UR, 2nd-order robust)",
                 "6H (0° UR, 3rd-order robust)"):
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


def test_main_window_facts_plots_and_stability():
    """Analyzer main window: no vector length field, InpoFact 1, curve data with the editor's facts, the
    Pulse_Sequence tab with text and plots, and the Stability tab = the fidelity panel."""
    from psa.GUI_pyside6 import PulseSequenceAnalyzerQt
    from psa.design import metrics
    from psa.design.fidelity import infidelity_grid
    main = PulseSequenceAnalyzerQt()
    main.show(); process(20)
    assert not hasattr(main, "vector_length_spin") and main.inpo_fact_spin.value() == 1
    path = os.path.join(tempfile.mkdtemp(), "p.csv")
    ps = PulseSequence(np.r_[np.full(60, 100.0), np.full(40, 50.0)], np.r_[np.zeros(60), np.full(40, 90.0)], 0.5e-6, 1e4)
    ps.save(path)
    main.pulse_sequence_edit.setText(path); main.update_analysis()
    assert main.latest_result is not None, main.info_error_text.toPlainText()
    text = main.curve_data_text.toPlainText()
    g = metrics.gate(ps)
    assert "100 sub-pulses × 0.5 µs = 50.00 µs" in text and f"{g['angle_deg']:.4f}°" in text
    assert "avg. power: 70.0 % of a rect. pulse" in text and "|A1|/T" in text and "vector area" in text
    cur = main.current_pulse()
    assert cur.N == 100 and np.allclose(cur.amplitude, ps.amplitude) and cur.umax == 1e4
    assert main.ps_matrix_text.toPlainText().splitlines()[2].split()[2:] == ["100.000000", "0.000000"]
    assert main.tabs.currentWidget() is main.error_plot and main.error_plot not in main._pending_plots
    assert main.trajectory_plot in main._pending_plots and main.torsion_plot in main._pending_plots
    main.tabs.setCurrentWidget(main.trajectory_plot)                # drawn when shown
    assert main.trajectory_plot not in main._pending_plots and len(main.trajectory_plot.ax.lines) == 1
    main._draw_curve_plots(all_plots=True); assert not main._pending_plots
    assert main._plots_dirty                                       # plots only drawn when visible
    main.tabs.setCurrentWidget(main.ps_subtabs); main.ps_subtabs.setCurrentWidget(main.ps_plot)
    assert len(main.ps_plot.figure.axes) == 4 and not main._plots_dirty
    main.offset_spin.setValue(1.0); main.update_analysis()          # same pulse: no redraw
    assert not main._plots_dirty
    main.offset_spin.setValue(0.0); main.tabs.setCurrentIndex(0)
    sp = main.stability_panel
    from psa.editor_gui.fidelity_window import AUTO, GOAL, OWN
    assert sp.ref_box.currentIndex() == OWN and sp.reference()[0] == {"kind": "own"}
    assert all(sp.ref_box.view().isRowHidden(i) and not sp.ref_box.model().item(i).isEnabled() for i in (AUTO, GOAL))
    assert not sp.free_widget.isHidden()                                   # freedom row of the own reference
    sp.map_off[2].setValue(11); sp.map_amp[2].setValue(3); sp.cut_off_range[2].setValue(21); sp.cut_amp_range[2].setValue(5)
    sp.refresh_button.click(); wait_until(lambda: not sp.busy, 120)
    assert sp.result is not None and sp.result["map"].shape == (3, 11), sp.status.text()
    assert np.allclose(sp.result["map"], infidelity_grid(ps, np.linspace(-10e3, 10e3, 11), [0.8, 1.0, 1.2]))
    main.offset_spin.setValue(2.0); main.update_analysis()
    assert "Shown curve (2 kHz, 100 %)" in main.curve_data_text.toPlainText()
    assert main.current_pulse() is cur and "has changed" not in sp.status.text()       # same pulse
    main.time_per_pulse_spin.setValue(1.0); main.update_analysis()
    assert main.current_pulse() is not cur and "the pulse has changed since" in sp.status.text()
    proc = sp._proc
    main.close(); process(20)
    assert sp._proc is None and not proc.is_alive()


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


def test_grape_tab_optimizes_robust_inversion():
    """GRAPE tab: rectangular start, local first-order robustness + 180°x gate, runs in a process."""
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.source_box.setCurrentText("rectangular")
    tab.dur_spin.setValue(150.0); tab.amp_spin.setValue(100.0); tab.phase_spin.setValue(0.0)
    tab.load_button.click(); process(20)
    assert tab.start_pulse.N == 300 and w.panel.pulse is tab.start_pulse
    cfg = tab.config()
    assert cfg.mode == "local" and cfg.orders == [1] and cfg.gate == "UR" and cfg.offsets_hz == [0.0]
    tab.solve_button.click()
    assert tab.solving and not tab.gate_box.isEnabled()
    wait_until(lambda: not tab.solving, 120)
    assert tab.status.text().startswith("all goals reached"), tab.status.text()
    out = tab.current_pulse()
    assert metrics.robustness(out)["A1_norm"] < 1e-8 and np.max(out.amplitude) <= 100 + 1e-9
    assert abs(metrics.gate(out)["angle_deg"] - 180) < 1e-4
    wait_panel(w)
    assert w.panel.pulse is out and tab.gate_box.isEnabled() and tab.solve_button.text() == "Optimize"
    tab.edit_button.click(); process(20)                         # hand over to the waveform tab
    assert w.current_tab() is w.waveform_tab and w.waveform_tab.current_pulse().N == out.N
    w.close()


def test_grape_tab_inputs_start_sources_and_stop():
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    # start from another tab
    w.text_tab.editor.setPlainText("R(50, 0, 0.5); R(50, 90, 0.5); R(50, 180, 0.5);"); w.text_tab.parse()
    tab.source_box.setCurrentText("Text tab"); tab.load_button.click(); process(20)
    assert tab.start_pulse.N == 3 and np.allclose(tab.start_pulse.phase, [0, 90, 180])
    # invalid inputs give readable errors instead of starting
    tab.offsets_edit.setText("1, x"); tab.solve_button.click()
    assert not tab.solving and "offsets" in tab.status.text()
    tab.offsets_edit.setText("-5:5:11")
    tab.order_checks[1].setChecked(False); tab.gate_box.setCurrentIndex(0); tab.solve_button.click()
    assert not tab.solving and "goal" in tab.status.text()
    # gate only = broadband PP config on the offset x B1 grid of the fields
    tab.gate_box.setCurrentIndex(2)
    assert tab.pp_widget.isVisibleTo(tab) and not tab.ur_widget.isVisibleTo(tab)
    tab.b1_edit.setText("90:110:3")
    assert tab.points_info.text().startswith("11 offsets × 3 B1 values = 33 points")
    tab.reg_box.setCurrentIndex(1); tab.reg_weight.setText("0.01")
    cfg = tab.config()
    assert cfg.mode == "broadband" and cfg.gate == "PP" and len(cfg.ensemble()) == 33
    assert np.allclose(cfg.ensemble()[:3], [(-5e3, 0.9), (-5e3, 1.0), (-5e3, 1.1)])
    assert cfg.regularization == "stay_close" and cfg.regularization_weight == 0.01
    tab.reweight_spin.setValue(4); tab.reweight_power_spin.setValue(0.8)
    cfg = tab.config()
    assert cfg.reweight_rounds == 4 and cfg.reweight_power == 0.8 and cfg.n_points() == 33
    tab.reweight_spin.setValue(0)
    # long random start: Stop returns the best pulse so far, the GUI stays usable
    tab.source_box.setCurrentText("random"); tab.dur_spin.setValue(400.0); tab.load_button.click(); process(20)
    w.solver_settings.stop_when_solved = False
    tab.solve_button.click()
    wait_until(lambda: tab._history or not tab.solving, 60)
    tab.solve_button.click()                                    # Stop
    wait_until(lambda: not tab.solving, 30)
    assert tab.status.text().startswith(("stopped", "all goals", "converged")), tab.status.text()
    assert tab.result is not None and tab.result.N == 800 and np.max(tab.result.amplitude) <= 100 + 1e-9
    before = tab.result
    tab.restart_button.click()
    assert tab.start_pulse.N == 800 and np.allclose(tab.start_pulse.amplitude, before.amplitude) and tab.result is None
    w.close()


def test_grape_tab_multistart_spline_and_limits():
    """K starts run in parallel processes (one residual line each); the best result is kept."""
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.source_box.setCurrentText("rectangular"); tab.dur_spin.setValue(150.0); tab.load_button.click(); process(20)
    tab.spline_radio.setChecked(True); tab.spline_spin.setValue(400)
    tab.solve_button.click()
    assert not tab.solving and "control points" in tab.status.text()             # more points than sub-pulses
    tab.spline_spin.setValue(16)
    tab.smooth_edit.setText("1e-5"); tab.amp_rate_spin.setValue(30); tab.phase_rate_spin.setValue(20)
    tab.starts_spin.setValue(3)
    cfg = tab.config()
    assert (cfg.controls, cfg.n_spline, cfg.smooth_weight, cfg.max_amp_rate, cfg.max_phase_rate_khz, cfg.n_starts) == \
        ("spline", 16, 1e-5, 30.0, 20.0, 3)
    tab.solve_button.click()
    assert tab.solving and not tab.spline_spin.isEnabled()
    wait_until(lambda: not tab.solving, 180)
    text = tab.status.text()
    assert text.startswith("all goals reached") and "best of 3 candidates" in text, text
    assert len(tab.conv_curves) == 3 and len(tab.last_infos) == 3
    out = tab.current_pulse()
    assert metrics.robustness(out)["A1_norm"] < 1e-8 and abs(metrics.gate(out)["angle_deg"] - 180) < 1e-4
    assert tab.spline_spin.isEnabled()
    w.close()


def test_grape_tab_exploration_presets_and_search():
    """Presets set the exploration inputs (correlation scales with the duration), editing one gives
    'custom'; a run with screening rounds and kicks reports the search."""
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.source_box.setCurrentText("rectangular"); tab.dur_spin.setValue(150.0); tab.load_button.click(); process(20)
    assert tab.preset_box.currentText() == "exploit" and tab.config().n_starts == 1
    tab.preset_box.setCurrentText("explore")
    cfg = tab.config()
    assert (cfg.n_starts, cfg.start_pool, cfg.screen_rounds, cfg.n_kicks, cfg.coarse_levels) == (24, "mix", 3, 8, [8, 16])
    assert abs(cfg.noise_correlation_us - 15.0) < 1e-9
    tab.advanced_button.click(); process(20)
    assert tab.explore_dialog.isVisible()
    tab.explore_dialog.kicks_spin.setValue(2)
    assert tab.preset_box.currentText() == "custom" and tab.config().n_kicks == 2
    tab.explore_dialog.coarse_edit.setText("8, x")
    tab.solve_button.click()
    assert not tab.solving and "coarse" in tab.status.text()
    tab.explore_dialog.hide()
    # a small search: 4 starts, 1 screening round (4 -> 2), 2 kicks (smoothing on, so kicks are useful)
    tab.preset_box.setCurrentText("balanced")
    tab.starts_spin.setValue(4); tab.explore_dialog.rounds_spin.setValue(1); tab.explore_dialog.screen_spin.setValue(8)
    tab.explore_dialog.kicks_spin.setValue(2); tab.smooth_edit.setText("1e-4")
    tab.solve_button.click()
    assert tab.solving and not tab.preset_box.isEnabled() and not tab.explore_dialog.kicks_spin.isEnabled()
    wait_until(lambda: not tab.solving, 240)
    text = tab.status.text()
    assert text.startswith("all goals reached") and "best of 6 candidates" in text and "round 1: 4 candidates" in text, text
    assert "kicks 1–2" in text and len(tab.last_infos) == 6 and tab.preset_box.isEnabled()
    out = tab.current_pulse()
    assert metrics.robustness(out)["A1_norm"] < 1e-8 and abs(metrics.gate(out)["angle_deg"] - 180) < 1e-4
    w.close()


def test_grape_tab_points_and_goals():
    """One point selection (offset × B1 grid of the fields, or custom points with weights); gate only maps to
    the broadband engine, gate + orders to the local one."""
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.offsets_edit.setText("-3, 0, 3"); tab.b1_edit.setText("95, 105")
    assert tab.points_info.text().startswith("3 offsets × 2 B1 values = 6 points")
    cfg = tab.config()
    assert cfg.mode == "local" and cfg.orders == [1] and cfg.n_points() == 6
    assert cfg.point_list()[1] == (-3e3, 1.05, 1.0)                 # local orders now also at B1 ≠ 100 %
    tab.order_checks[1].setChecked(False)
    assert tab.config().mode == "broadband" and tab.config().n_points() == 6
    tab.b1_edit.setText("95, 0")
    assert "positive" in tab.points_info.text()
    tab.solve_button.click()
    assert not tab.solving and "positive" in tab.status.text()
    # custom selection with weights: fields locked, back to grid restores them
    tab.b1_edit.setText("100")
    tab.set_points([[-2e3, 1.0, 1.0], [1e3, 0.9, 4.0]])
    assert not tab.offsets_edit.isEnabled() and tab.grid_button.isVisibleTo(tab)
    assert tab.points_info.text() == "custom selection: 2 points, weights 1 … 4"
    cfg = tab.config()
    assert cfg.point_list() == [(-2e3, 1.0, 1.0), (1e3, 0.9, 4.0)] and np.allclose(cfg.user_weights(), [0.4, 1.6])
    tab.grid_button.click()
    assert tab.custom_points is None and tab.offsets_edit.isEnabled() and tab.offsets_edit.text() == "-3, 0, 3"
    # a selection that is a full grid with weights 1 goes back into the fields
    tab.set_points([[o * 1e3, b, 1.0] for o in np.linspace(-10, 10, 21) for b in (0.9, 1.0)])
    assert tab.custom_points is None and tab.offsets_edit.text() == "-10:10:21" and tab.b1_edit.text() == "90, 100"
    # a short gate-only run reports the score over the area (between the points too)
    tab.source_box.setCurrentText("rectangular"); tab.dur_spin.setValue(150.0); tab.load_button.click(); process(20)
    tab.offsets_edit.setText("-1, 1"); tab.b1_edit.setText("100"); tab.angle_spin.setValue(90.0)
    assert tab.area_check.isChecked()
    tab.solve_button.click()
    wait_until(lambda: not tab.solving, 120)
    text = tab.status.text()
    assert "area -1 … 1 kHz (201×1 test points): worst infidelity" in text, text
    assert tab.last_area is not None and tab.last_area["grid"] == (1, 201)
    w.close()


def test_grape_point_editor():
    """Visual point editor: add (mouse and API), mirror, lock to grid, weights, fill, select / move / delete,
    undo / redo, save / load, background map, Apply to the tab."""
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.offsets_edit.setText("-2, 2"); tab.b1_edit.setText("100")
    tab.visual_button.click(); process(50)
    ed = tab.point_editor
    assert ed.isVisible() and ed.points() == [[-2e3, 1.0, 1.0, None], [2e3, 1.0, 1.0, None]]
    # a real click on empty space adds a point there
    ed.view.setRange(xRange=(-10, 10), yRange=(80, 120), padding=0); process(50)
    QTest.mouseClick(ed.plot.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     view_to_widget(ed.plot, 5.0, 110.0)); process(20)
    assert len(ed.pts) == 3 and abs(ed.pts[2][0] - 5.0) < 0.3 and abs(ed.pts[2][1] - 110) < 1.5 and ed.sel == {2}
    # lock to grid snaps (also the existing points), mirror adds −offset
    w.solver_settings.point_grid_offset_khz, w.solver_settings.point_grid_b1_percent = 0.5, 1.0
    ed.lock_check.setChecked(True)
    assert all(abs(p[0] * 2 - round(p[0] * 2)) < 1e-9 and abs(p[1] - round(p[1])) < 1e-9 for p in ed.pts)
    ed.mirror_check.setChecked(True); ed.weight_spin.setValue(3.0)
    ed.add_point(3.14, 95.3)
    assert [3.0, 95.0, 3.0, None] in ed.pts and [-3.0, 95.0, 3.0, None] in ed.pts and len(ed.sel) == 2
    n = len(ed.pts)
    ed.undo(); assert len(ed.pts) == n - 2
    ed.redo(); assert len(ed.pts) == n
    # select in a rectangle, move (snapped), weight, delete
    ed.select_in_rect(-4, 90, 4, 99)
    assert ed.sel == {i for i, p in enumerate(ed.pts) if p[1] == 95.0}
    ed.move_selected(0.26, 2.2)
    assert sorted(ed.pts[i][1] for i in ed.sel) == [97.0, 97.0] and {ed.pts[i][0] for i in ed.sel} == {-2.5, 3.5}
    ed.set_weight_selected(0.5)
    assert all(ed.pts[i][2] == 0.5 for i in ed.sel)
    ed.delete_selected(); assert len(ed.pts) == n - 2
    # fill a rectangle with a grid (no mirror)
    ed.mirror_check.setChecked(False); ed.weight_spin.setValue(1.0)
    ed.fill_rect(-1, 90, 1, 110, n_off=3, n_b1=2)
    assert {(p[0], p[1]) for p in ed.pts} >= {(-1.0, 90.0), (0.0, 110.0), (1.0, 110.0)}
    # table editing
    r = ed.pts.index([-1.0, 90.0, 1.0, None])
    ed.table.item(r, 2).setText("7"); process(5)
    assert ed.pts[r][2] == 7.0
    # save / load round trip (JSON and CSV)
    with tempfile.TemporaryDirectory() as d:
        before = ed.points()
        for name in ("sel.json", "sel.csv"):
            ed.save(os.path.join(d, name)); ed.delete_selected(); ed.select_all(); ed.delete_selected()
            assert not ed.pts and not ed.apply_button.isEnabled()
            ed.load(os.path.join(d, name))
            assert np.allclose([p[:3] for p in ed.points()], [p[:3] for p in before]), name
            assert [p[3] for p in ed.points()] == [p[3] for p in before], name
    # background: infidelity of the current pulse for the gate goal
    offs, b1, z = ed.background_map(21, 5)
    assert z.shape == (5, 21) and np.all(np.isfinite(z)) and np.all(z <= 1e-9)
    ed.bg_check.setChecked(True)
    assert ed.bg_item.isVisible()
    tab.gate_box.setCurrentIndex(0)
    assert ed.background_map() is None
    tab.gate_box.setCurrentIndex(1)
    # Apply: custom selection in the tab, used by the config
    ed.apply_button.click(); process(10)
    assert tab.custom_points is not None and not tab.offsets_edit.isEnabled()
    cfg = tab.config()
    assert cfg.n_points() == len(ed.pts) and cfg.user_weights() is not None
    assert sorted(cfg.point_list()) == sorted((p[0] * 1e3, p[1] / 100, p[2]) for p in ed.pts)
    ed.close(); w.close()


def test_fidelity_window():
    """⊕ on the panel opens the fidelity window; Refresh computes map and cuts in a worker process (a new
    Refresh replaces a running one), also while a GRAPE optimization runs; pulse changes mark it outdated."""
    from psa.design.fidelity import infidelity_grid
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.source_box.setCurrentText("rectangular"); tab.dur_spin.setValue(150.0); tab.load_button.click(); process(20)
    wait_panel(w)
    w.panel.fid_expand.click(); process(50)
    fw = w.panel.fidelity_window
    assert fw.isVisible() and fw.result is None and "Refresh" in fw.status.text()
    fw.map_off[2].setValue(21); fw.map_amp[2].setValue(7); fw.cut_off_range[2].setValue(51); fw.cut_amp_range[2].setValue(31)
    # invalid range: readable error, nothing starts
    fw.map_amp[0].setValue(130.0); fw.refresh_button.click()
    assert not fw.busy and "range" in fw.status.text()
    fw.map_amp[0].setValue(80.0)
    # a big job replaced at once by a small one: only the second result arrives
    fw.map_off[2].setValue(2001); fw.refresh_button.click(); assert fw.busy
    fw.map_off[2].setValue(21); fw.refresh_button.click()
    wait_until(lambda: not fw.busy, 120)
    r = fw.result
    assert r is not None and r["id"] == 2 and r["map"].shape == (7, 21), fw.status.text()
    assert r["cut_offset"].shape == (51,) and r["cut_amp"].shape == (31,) and fw.status.text().startswith("computed")
    ps = tab.start_pulse
    assert np.allclose(r["map"], infidelity_grid(ps, np.linspace(-10e3, 10e3, 21), np.linspace(0.8, 1.2, 7)))
    x, y = fw.cut_off_curve.getData()
    assert len(x) == 51 and np.allclose(y, np.log10(np.maximum(r["cut_offset"], 1e-16)))
    fw.quantity_box.setCurrentIndex(1)                     # F instead of log10(1 - F): no recompute
    assert np.allclose(fw.cut_off_curve.getData()[1], 1 - r["cut_offset"])
    fw.view_box.setCurrentIndex(1); process(20)            # heat map
    assert fw.map_info.text().startswith("worst 1 − F")
    # references: automatic = the GRAPE gate goal while the GRAPE tab is shown / optimizing, else the own gate
    from psa.editor_gui.fidelity_window import AUTO, GOAL, OWN, PP, UR
    assert fw.ref_box.currentIndex() == AUTO
    tab.angle_spin.setValue(90.0); tab.axis_spins[0].setValue(0.0); tab.axis_spins[1].setValue(1.0)
    assert fw.reference()[0] == {"kind": "UR", "axis": [0.0, 1.0, 0.0], "angle": 90.0}
    assert "GRAPE gate goal: UR 90° about (0, 1, 0)" in fw.reference()[1]
    w.show_tab(w.text_tab)
    assert fw.reference()[0] == {"kind": "own"}
    fw.ref_box.setCurrentIndex(GOAL)                       # explicit: the GRAPE goal from any tab
    assert fw.reference()[0]["angle"] == 90.0
    tab.gate_box.setCurrentIndex(2)                        # PP goal
    assert fw.reference()[0] == {"kind": "PP", "m0": [0.0, 0.0, 1.0], "mt": [1.0, 0.0, 0.0]}
    tab.gate_box.setCurrentIndex(0)
    try:
        fw.reference()
    except ValueError:
        pass
    else:
        raise AssertionError("GRAPE goal without a gate goal accepted")
    fw.refresh_button.click()
    assert not fw.busy and "no gate goal" in fw.status.text()
    tab.gate_box.setCurrentIndex(1)
    fw.ref_box.setCurrentIndex(UR)                         # the UR inputs start from the GRAPE goal
    assert fw.angle_spin.value() == 90.0 and fw.axis_spins[1].value() == 1.0
    fw.angle_spin.setValue(45.0)
    assert fw.reference()[0]["angle"] == 45.0 and fw.reference()[1] == "UR 45° about (0, 1, 0)"
    fw.ref_box.setCurrentIndex(PP)
    assert fw.reference()[0]["kind"] == "PP"
    fw.ref_box.setCurrentIndex(OWN)
    assert fw.reference()[0] == {"kind": "own"}
    fw.ref_box.setCurrentIndex(AUTO)
    w.show_tab(tab)
    # refresh while GRAPE optimizes: both finish; automatic uses the running optimization's gate goal
    w.solver_settings.stop_when_solved = False
    tab.solve_button.click(); assert tab.solving
    w.show_tab(w.text_tab)                                 # still the GRAPE goal: it is optimizing
    fw.refresh_button.click()
    wait_until(lambda: not fw.busy, 120)
    assert fw.result["id"] == 3 and fw.result["job"]["reference"] == {"kind": "UR", "axis": [0.0, 1.0, 0.0], "angle": 90.0}
    assert "reference: GRAPE gate goal: UR 90°" in fw.status.text() and tab.solving
    w.show_tab(tab)                                        # the panel follows the shown tab's pulse
    wait_until(lambda: tab._history or not tab.solving, 60)
    tab.solve_button.click(); wait_until(lambda: not tab.solving, 30)
    wait_panel(w)
    assert "the pulse has changed since" in fw.status.text()
    proc = fw._proc
    fw.close(); process(20)
    assert fw._proc is None and not proc.is_alive()
    w.close()


def test_fidelity_window_freedom():
    """Freedom of the window's own UR / PP reference: fitted on the map area or fixed; the values are shown."""
    from psa.design.fidelity import infidelity_grid
    from psa.editor_gui.fidelity_window import AUTO, GOAL, PP, UR
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.source_box.setCurrentText("rectangular"); tab.dur_spin.setValue(100.0); tab.load_button.click(); process(20)
    wait_panel(w)
    w.panel.fid_expand.click(); process(50)
    fw = w.panel.fidelity_window
    fw.map_off[2].setValue(21); fw.map_amp[2].setValue(5); fw.cut_off_range[2].setValue(41); fw.cut_amp_range[2].setValue(11)
    for i in (AUTO, GOAL):                                 # the GRAPE tab's freedom there: no own row
        fw.ref_box.setCurrentIndex(i); assert fw.free_widget.isHidden()
    fw.ref_box.setCurrentIndex(PP)
    fw.set_reference_inputs({"kind": "PP", "m0": [0, 0, 1], "mt": [1, 0, 0]})
    assert not fw.free_widget.isHidden() and fw.free_source.isHidden()
    fw.free_box.setCurrentIndex(2)                         # linear phase, fitted: max |delay| shown
    assert not fw.max_delay_spin.isHidden() and fw.tau_spin.isHidden() and fw.phi_spin.isHidden()
    ref, label = fw.reference()
    assert ref == {"kind": "PP", "m0": [0.0, 0.0, 1.0], "mt": [1.0, 0.0, 0.0], "freedom": "linear", "fit": True,
                   "max_delay_us": 0.0}
    assert "up to a linear phase (phase + delay) (fitted on the map area)" in label
    fw.refresh_button.click(); wait_until(lambda: not fw.busy, 120)
    c = fw.result["correction"]
    assert c["source"] == "fitted on the map area" and "correction: phase" in fw.status.text() and "delay" in fw.status.text()
    plain = infidelity_grid(tab.start_pulse, np.linspace(-10e3, 10e3, 21), np.linspace(0.8, 1.2, 5),
                            {"kind": "PP", "m0": [0, 0, 1], "mt": [1, 0, 0]})
    assert np.mean(fw.result["map"]) <= np.mean(plain) + 1e-15
    fitted = fw.result["map"]
    fw.set_correction(c["phi_deg"], c["tau_us"])           # fixed = the fitted values: the same map
    assert not fw.tau_spin.isHidden() and fw.max_delay_spin.isHidden()
    assert "fixed:" in fw.reference()[1]
    fw.refresh_button.click(); wait_until(lambda: not fw.busy, 120)
    assert fw.result["correction"]["source"] == "fixed" and np.allclose(fw.result["map"], fitted, atol=1e-3)   # spins round
    fw.map_check.setChecked(False); fw.free_source.setCurrentIndex(0)   # no map: still fitted on the map ranges
    fw.refresh_button.click(); wait_until(lambda: not fw.busy, 120)
    assert fw.result["map"] is None and np.isclose(fw.result["correction"]["tau_us"], c["tau_us"])
    fw.free_box.setCurrentIndex(0)
    assert "freedom" not in fw.reference()[0]
    fw.close(); w.close()


def test_grape_quick_b1_range():
    """Quick mode B1 range: 100 to 100 % = off (plain recipe); a range adds the B1 values of the rule."""
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.quick_radio.setChecked(True)
    tab.q_dur_spin.setValue(450.0); tab.q_lo_spin.setValue(-10.0); tab.q_hi_spin.setValue(10.0)
    assert tab.quick_b1_range() == (1.0, 1.0) and tab.q_points_info.text() == "19 offsets (nominal amplitude)"
    assert tab.quick_setup()[1].bb_scale_n == 1
    tab.q_b1_lo_spin.setValue(80.0); tab.q_b1_hi_spin.setValue(120.0)
    assert "19 offsets × 5 B1 values = 95 points" in tab.q_points_info.text()
    cfg = tab.quick_setup()[1]
    assert (cfg.bb_scale_min, cfg.bb_scale_max, cfg.bb_scale_n) == (0.8, 1.2, 5)
    tab.q_b1_hi_spin.setValue(70.0)
    try:
        tab.quick_setup()
    except ValueError as e:
        assert "B1" in str(e)
    else:
        raise AssertionError("accepted an inverted B1 range")
    w.close()


def test_grape_per_point_goals():
    """Selective pulses: own goals per point in the point editor (default: the GRAPE gate goal), shown in colours
    and the table, used by the config, the optimization, the area score and the fidelity window."""
    from psa.editor_gui.goal_dialog import PRESETS, GoalDialog
    ex = {"gate": "PP", "m0": [0, 0, 1], "mt": [1, 0, 0]}
    keep = {"gate": "PP", "m0": [0, 0, 1], "mt": [0, 0, 1]}
    # the goal dialog: presets, default, validation
    dlg = GoalDialog(None, 3)
    assert dlg.goal() is None and dlg.default_radio.isChecked()
    dlg._preset([p[0] for p in PRESETS].index("excite z → x (PP)"))
    assert dlg.goal() == {"gate": "PP", "m0": [0.0, 0.0, 1.0], "mt": [1.0, 0.0, 0.0]}
    dlg.ur_radio.setChecked(True); dlg.angle_spin.setValue(30.0)
    assert dlg.goal()["angle"] == 30.0
    for sp in dlg.axis_spins:
        sp.setValue(0.0)
    dlg._ok()
    assert dlg.result() != 1 and "non-zero" in dlg.error.text()
    dlg.close()
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.source_box.setCurrentText("rectangular"); tab.dur_spin.setValue(300.0); tab.amp_spin.setValue(30.0)
    tab.load_button.click(); process(20)
    tab.order_checks[1].setChecked(False); tab.gate_box.setCurrentIndex(2)          # PP z -> x default
    tab.offsets_edit.setText("-1:1:3, -7:-5:3, 5:7:3"); tab.b1_edit.setText("100")
    tab.visual_button.click(); process(50)
    ed = tab.point_editor
    assert all(p[3] is None for p in ed.pts) and "default: PP z→x" in ed.legend.text()
    tab.mt_spins[0].setValue(0.0); tab.mt_spins[1].setValue(1.0)                    # legend follows the gate goal
    assert "default: PP z→y" in ed.legend.text()
    tab.mt_spins[1].setValue(0.0); tab.mt_spins[0].setValue(1.0)
    ed.select_in_rect(4, 90, 8, 110); ed.select_in_rect(-8, 90, -4, 110, extend=True)
    assert len(ed.sel) == 6
    ed.set_goal_selected(keep)
    assert sum(p[3] == keep for p in ed.pts) == 6 and "PP z→z" in ed.legend.text() and "default" in ed.legend.text()
    assert {ed.table.item(r, 3).text() for r in range(len(ed.pts))} == {"default", "PP z→z"}
    assert "6 with own goal" in ed.info.text()
    ed.undo(); assert all(p[3] is None for p in ed.pts)
    ed.redo(); assert sum(p[3] is not None for p in ed.pts) == 6
    offs, b1, z = ed.background_map(31, 3)                 # nearest point's goal
    assert np.all(np.isfinite(z))
    ed.apply_button.click(); process(10)
    assert tab.custom_points is not None and "6 with own goal" in tab.points_info.text()
    cfg = tab.config()
    assert cfg.mode == "broadband" and sum(g is not None for g in cfg.point_goals()) == 6
    by_offset = dict(zip([o for o, _ in cfg.ensemble()], cfg.point_targets()))
    assert by_offset[0.0]["mt"] == [1.0, 0.0, 0.0] and by_offset[7e3] == keep and by_offset[-5e3] == keep
    # gate "none" works too when every point has its own goal
    ed.select_all(); ed.set_goal_selected(ex); ed.select_in_rect(4, 90, 8, 110); ed.set_goal_selected(keep)
    ed.apply_button.click(); process(10)
    tab.gate_box.setCurrentIndex(0)
    assert tab.config().gate == "none"
    tab.gate_box.setCurrentIndex(2)
    # a short selective optimization through the tab; area score and fidelity window use the goals per point
    w.solver_settings.max_nfev = 150
    tab.solve_button.click()
    wait_until(lambda: not tab.solving, 180)
    assert tab.last_area is not None and "area -7 … 7 kHz" in tab.status.text(), tab.status.text()
    wait_panel(w)
    w.panel.fid_expand.click(); process(50)
    fw = w.panel.fidelity_window
    ref, label = fw.reference()
    assert ref["kind"] == "points" and "goals per point (9 points, 9 with own goal" in label
    fw.map_off[2].setValue(15); fw.map_amp[2].setValue(3); fw.cut_off_range[2].setValue(31); fw.cut_amp_range[2].setValue(5)
    fw.refresh_button.click(); wait_until(lambda: not fw.busy, 120)
    assert fw.result is not None and fw.result["map"].shape == (3, 15) and np.all(np.isfinite(fw.result["map"]))
    ed.close(); fw.close(); w.close()


def test_grape_freedom():
    """Goals up to a free phase / linear phase: settings in the gate group (hidden in quick mode), the fitted
    values in the status, the fidelity window judges with the same fitted rotation."""
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    assert tab.freedom == "none" and not tab.delay_spin.isVisibleTo(tab)
    tab.freedom_box.setCurrentIndex(2)
    assert tab.freedom == "linear" and tab.delay_spin.isVisibleTo(tab)
    tab.quick_radio.setChecked(True)
    assert not tab.freedom_widget.isVisibleTo(tab) and tab.quick_setup()[1].freedom == "none"
    tab.expert_radio.setChecked(True)
    assert tab.freedom_widget.isVisibleTo(tab)
    tab.source_box.setCurrentText("rectangular"); tab.dur_spin.setValue(150.0); tab.load_button.click(); process(20)
    tab.order_checks[1].setChecked(False); tab.gate_box.setCurrentIndex(2)          # PP z -> x
    tab.offsets_edit.setText("-8:8:17"); tab.b1_edit.setText("100")
    cfg = tab.config()
    assert cfg.freedom == "linear" and cfg.max_delay_us == 0.0 and cfg.mode == "broadband"
    w.solver_settings.max_nfev = 150
    tab.solve_button.click(); assert not tab.freedom_box.isEnabled()
    wait_until(lambda: not tab.solving, 180)
    text = tab.status.text()
    assert "free phase [°]" in text and "free delay [µs]" in text and tab.freedom_box.isEnabled(), text
    wait_panel(w)
    w.panel.fid_expand.click(); process(50)
    fw = w.panel.fidelity_window
    ref, label = fw.reference()
    assert ref["kind"] == "points" and ref["freedom"] == "linear" and "up to a linear phase" in label
    fw.map_check.setChecked(False); fw.cut_amp_check.setChecked(False)
    fw.cut_off_range[0].setValue(-8.0); fw.cut_off_range[1].setValue(8.0); fw.cut_off_range[2].setValue(17)
    fw.refresh_button.click(); wait_until(lambda: not fw.busy, 120)
    from psa.design.grape import GrapeProblem
    prob = GrapeProblem(tab.result, tab._cfg)
    summ = prob.goals(prob.x0)[2]
    assert abs(np.max(fw.result["cut_offset"]) - summ["max infidelity"]) < 1e-9 * max(1.0, summ["max infidelity"])
    fw.close(); w.close()


def test_grape_tab_quick_mode():
    """Quick mode: only gate, duration, band, limits and time budget; runs the study recipe for the budget and
    keeps the best pulse within the rate limits."""
    from psa.design.grape import band_quality, pulse_rates
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    assert not tab.quick_box.isVisibleTo(tab) and tab.expert_boxes[0].isVisibleTo(tab)
    tab.gate_box.setCurrentIndex(0)
    tab.quick_radio.setChecked(True)
    assert tab.quick_box.isVisibleTo(tab) and not any(b.isVisibleTo(tab) for b in tab.expert_boxes)
    assert tab.gate_box.currentIndex() == 1 and tab.gate_box.isVisibleTo(tab)         # quick mode needs a gate
    tab.q_lo_spin.setValue(5.0); tab.q_hi_spin.setValue(-5.0)
    tab.solve_button.click()
    assert not tab.solving and "bandwidth" in tab.status.text()
    tab.q_lo_spin.setValue(-4.0); tab.q_hi_spin.setValue(4.0)
    tab.q_dur_spin.setValue(100.0); tab.angle_spin.setValue(90.0); tab.q_budget_spin.setValue(0.5)
    ps0, cfg, st = tab.quick_setup()
    assert ps0.N == 200 and cfg.bb_offset_n == 3 and cfg.n_starts == 22 and st.max_nfev == 3000
    tab.solve_button.click()
    assert tab.solving and not tab.q_dur_spin.isEnabled() and not tab.expert_radio.isEnabled()
    wait_until(lambda: not tab.solving, 120)
    text = tab.status.text()
    assert text.startswith("quick run finished") and "worst infidelity" in text, text
    assert tab.last_quick["within_limits"] and tab.last_quick["elapsed_s"] < 30 + 15
    out = tab.current_pulse()
    assert out.N == 200 and pulse_rates(out)[0] <= 21 and pulse_rates(out)[1] <= 52.5
    start_worst = band_quality(ps0, cfg, -4e3, 4e3)["worst infidelity"]
    assert band_quality(out, cfg, -4e3, 4e3)["worst infidelity"] < start_worst / 10
    wait_panel(w)
    assert tab.q_dur_spin.isEnabled() and w.panel.pulse is out
    tab.expert_radio.setChecked(True)
    assert tab.expert_boxes[0].isVisibleTo(tab) and not tab.quick_box.isVisibleTo(tab)
    w.close()


def test_grape_quick_series_library():
    """Quick mode output: a single pulse saved under the naming rule; a series (library) in its own folder with
    library.json, shown in the library window; a click loads a pulse."""
    import json
    import tempfile
    from psa.design.library import read_index
    d = tempfile.mkdtemp()
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.quick_radio.setChecked(True)
    tab.gate_box.setCurrentIndex(1); tab.angle_spin.setValue(90.0)
    for sp, v in zip(tab.axis_spins, (1.0, 0.0, 0.0)):
        sp.setValue(v)
    tab.q_lo_spin.setValue(-4.0); tab.q_hi_spin.setValue(4.0); tab.q_dur_spin.setValue(60.0)
    tab.q_budget_spin.setValue(0.5); tab.q_folder_edit.setText(d)
    assert tab.q_dur_hi_spin.isHidden() and tab.q_count_spin.isHidden()
    tab.solve_button.click(); wait_until(lambda: not tab.solving, 120)
    assert tab.last_saved == os.path.join(d, "pulse_0060us_UR100-90_off8_B10_1.bruker"), tab.status.text()
    assert os.path.exists(tab.last_saved) and "saved:" in tab.status.text()
    # series: 3 durations, total 0.5 min split in proportion; an existing name gets the next index
    tab.q_series_check.setChecked(True)
    assert not tab.q_dur_hi_spin.isHidden() and tab.q_budget_label.text() == "total time:"
    tab.q_dur_spin.setValue(40.0); tab.q_dur_hi_spin.setValue(60.0); tab.q_count_spin.setValue(3)
    assert not tab.q_warm_check.isHidden()
    tab.q_warm_check.setChecked(True)
    durs, budgets = tab.series_plan()
    assert [round(x * 1e6) for x in durs] == [40, 50, 60] and np.isclose(sum(budgets), 30.0)
    assert budgets[2] > budgets[0] and "3 pulses: 40, 50, 60 µs" in tab.q_points_info.text()
    os.makedirs(os.path.join(d, "pulse_0040-0060us_UR100-90_off8_B10_1"))
    tab.solve_button.click()
    lib = os.path.join(d, "pulse_0040-0060us_UR100-90_off8_B10_2")
    assert os.path.isdir(lib) and tab.library_window.isVisible()
    wait_until(lambda: not tab.solving and tab._series is None, 240)
    idx = read_index(lib)
    assert [p["file"] for p in idx["pulses"]] == ["pulse_0040us.bruker", "pulse_0050us.bruker", "pulse_0060us.bruker"]
    assert idx["settings"]["n_pulses"] == 3 and all(os.path.exists(os.path.join(lib, p["file"])) for p in idx["pulses"])
    from psa.design.quick import stretch_pulse                         # warm start: the 60 µs run began from 50 µs
    prev = PulseSequence.load(os.path.join(lib, "pulse_0050us.bruker"), tab.dt, tab.umax)
    assert idx["settings"]["warm_start"] and tab.start_pulse.N == 120
    assert np.allclose(tab.start_pulse.amplitude, stretch_pulse(prev, 60e-6).amplitude, atol=1e-4)
    assert tab.status.text().startswith("series finished: 3 of 3 pulses") and tab.last_library == lib
    lw = tab.library_window
    assert lw.table.rowCount() == 3 and len(lw.worst_curve.getData()[0]) == 3
    lw._row_clicked(0, 0)
    assert tab.result.N == 80                                           # 40 µs
    lw.set_library(None, [])
    assert lw.open_folder(lib) and lw.table.rowCount() == 3 and not lw.open_folder(d)
    tab.q_folder_edit.setText("")
    tab.solve_button.click()
    assert not tab.solving and "output folder" in tab.status.text()
    w.close()


def test_grape_stop_without_result_keeps_the_shown_pulse():
    """Stop when no candidate returned in time (busy machine): the last pulse shown is kept, not lost."""
    from psa.editor_gui.tab_grape import STOP_GRACE_S
    w = make_window()
    tab = w.grape_tab
    w.show_tab(tab)
    tab.source_box.setCurrentText("rectangular"); tab.dur_spin.setValue(50.0); tab.load_button.click(); process(20)
    shown = PulseSequence(np.full(100, 60.0), np.zeros(100), tab.dt, tab.umax)

    class NoResult:
        cands = {}
        def best(self):
            raise RuntimeError("no result")

    tab._quick, tab._search, tab._snapshot = None, NoResult(), shown
    tab._complete()
    assert tab.result is shown and tab.status.text().startswith("stopped") and f"{STOP_GRACE_S:g} s" in tab.status.text()
    assert STOP_GRACE_S == 10.0
    tab._snapshot = None
    tab._complete()
    assert "no result" in tab.status.text()
    w.close()


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
