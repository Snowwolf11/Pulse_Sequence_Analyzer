"""Tests of the pulse file names and library bookkeeping (psa.design.library)."""
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.grape import GrapeConfig
from psa.design.library import (config_goal, duration_text, goal_tag, int_vector, library_folder, number_text,
                                range_tag, read_index, series_budgets, series_durations, single_path, write_index)
from psa.design.quick import quick_config

DT, UMAX = 0.5e-6, 1e4


def test_numbers_round_half_up():
    assert duration_text(549.5e-6) == "0550us" and duration_text(405e-6) == "0405us" and duration_text(5e-6) == "0005us"
    assert duration_text(1234.4e-6) == "1234us"
    assert [number_text(x) for x in (20.0, 7.25, 0.05, 40.000000001, 2.5e-7)] == ["20", "7p3", "0p1", "40", "0"]
    assert number_text(2.5, 0) == "3" and number_text(-89.5, 0) == "m90" and number_text(-0.04) == "0"


def test_vectors_and_goals():
    assert int_vector([2, 3, 5]) == (2, 3, 5) and int_vector([0.2, 0.3, 0.5]) == (2, 3, 5)
    assert int_vector([4, 6, 10]) == (2, 3, 5) and int_vector([0, 0, 7]) == (0, 0, 1)
    assert int_vector([1, 1 / 3, 0]) == (3, 1, 0) and int_vector([0, -1, 0]) == (0, -1, 0)
    v = int_vector([0.1234, 0.9, 0.4])                                   # no exact small direction: the closest
    assert max(map(abs, v)) <= 9
    assert goal_tag({"gate": "UR", "axis": [2, 3, 5], "angle": 180}) == "UR235-180"
    assert goal_tag({"gate": "UR", "axis": [1, 0, 0], "angle": 89.5}) == "UR100-90"
    assert goal_tag({"gate": "PP", "m0": [0, 0, 1], "mt": [0, -1, 0]}) == "PP001-0m10"
    assert goal_tag({"gate": "PP", "m0": [0, 0, 1], "mt": [0.6, 0.8, 0]}) == "PP001-340"


def test_names_of_a_config():
    cfg = quick_config(550e-6, (-1e4, 1e4), "UR", (2, 3, 5), 180.0, dt=DT, b1_range=(0.8, 1.2), umax_hz=UMAX)[0]
    assert range_tag(cfg) == "off20_B140"
    # the goal most points have
    pp = {"gate": "PP", "m0": [0, 0, 1], "mt": [1, 0, 0]}
    sel = GrapeConfig(mode="broadband", gate="UR", target_axis=[1, 0, 0], target_angle_deg=90.0,
                      points=[[-2e3, 1.0, 1.0, pp], [0.0, 0.95, 1.0, pp], [2e3, 1.05, 1.0, None]])
    assert config_goal(sel) == pp and range_tag(sel) == "off4_B110"
    half = GrapeConfig(mode="broadband", gate="UR", points=[[-3.75e3, 0.975, 1.0], [3.75e3, 1.0, 1.0]])
    assert range_tag(half) == "off7p5_B12p5"
    with tempfile.TemporaryDirectory() as d:
        p1 = single_path(d, cfg, 549.6e-6)
        assert os.path.basename(p1) == "pulse_0550us_UR235-180_off20_B140_1.bruker"
        open(p1, "w").close()
        assert single_path(d, cfg, 550e-6).endswith("_2.bruker")              # never overwrite
        lib = library_folder(d, cfg, [100e-6, 600e-6])
        assert os.path.basename(lib) == "pulse_0100-0600us_UR235-180_off20_B140_1"
        os.makedirs(lib)
        write_index(lib, {"n_pulses": 2}, [{"duration_us": 100.0, "file": "pulse_0100us.bruker"}])
        assert read_index(lib)["pulses"][0]["file"] == "pulse_0100us.bruker"
        assert library_folder(d, cfg, [100e-6, 600e-6]).endswith("_2")


def test_series():
    durs = series_durations(100, 600, 11)
    assert [round(x * 1e6) for x in durs] == list(range(100, 601, 50))
    assert [round(x * 1e6) for x in series_durations(100, 101, 5)] == [100, 101]       # rounded, no duplicates
    b = series_budgets(durs, 600.0)
    assert np.isclose(sum(b), 600.0) and np.isclose(b[-1] / b[0], 6.0)


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t(); print(f"[PASS] {t.__name__}")
        except AssertionError as e:
            failed += 1; print(f"[FAIL] {t.__name__}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
