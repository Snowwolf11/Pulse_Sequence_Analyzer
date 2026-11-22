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
