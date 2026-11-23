"""Tests of the GRAPE point selections (psa.design.points)."""
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.points import format_values, grid_points, is_grid, load_points, parse_values, save_points, snap


def test_parse_values():
    assert parse_values("0") == [0.0]
    assert parse_values("-10:10:5, 12; 15 0") == [-10.0, -5.0, 0.0, 5.0, 10.0, 12.0, 15.0]   # duplicate 0 dropped
    assert parse_values("3:3:1") == [3.0]
    assert np.allclose(parse_values("90:110:3"), [90, 100, 110])
    for bad in ("", "1, x", "1:2", "1:2:0", "1:2:3:4", "inf"):
        try:
            parse_values(bad, "offsets")
        except ValueError as e:
            assert "offsets" in str(e)
            continue
        raise AssertionError(f"accepted {bad!r}")


def test_grid_and_format_round_trip():
    pts = grid_points([-2.0, 0.0, 2.0], [90.0, 100.0])
    assert pts[:2] == [[-2e3, 0.9, 1.0], [-2e3, 1.0, 1.0]] and len(pts) == 6
    assert is_grid(pts) == ([-2.0, 0.0, 2.0], [90.0, 100.0])
    assert is_grid(pts[:5]) is None                                  # not complete
    assert is_grid([[0.0, 1.0, 2.0]]) is None                        # weight != 1
    assert format_values(np.linspace(-10, 10, 21)) == "-10:10:21" and format_values([1, 2.5]) == "1, 2.5"
    assert parse_values(format_values(np.linspace(-7.5, 7.5, 31))) == list(np.linspace(-7.5, 7.5, 31))
    try:
        grid_points([0.0], [0.0])
    except ValueError:
        pass
    else:
        raise AssertionError("accepted B1 = 0 %")


def test_snap():
    assert snap(1.234, 97.26) == (1.2, 97.5)
    assert snap(-0.05001, 100.2, {"offset_khz": 0.5, "b1_percent": 1.0}) == (-0.0, 100.0)
    assert snap(0.123, 99.9, {"offset_khz": 0, "b1_percent": 0}) == (0.123, 99.9)          # 0 = no snapping


def test_save_and_load_json_and_csv():
    pts = [[-2500.0, 0.95, 1.0], [0.0, 1.0, 3.5], [4100.0, 1.1, 0.25]]
    with tempfile.TemporaryDirectory() as d:
        for name in ("sel.json", "sel.csv"):
            path = os.path.join(d, name)
            save_points(path, pts, {"offset_khz": 0.2, "b1_percent": 1.0})
            back, grid = load_points(path)
            assert np.allclose(back, pts), name
            assert (grid == {"offset_khz": 0.2, "b1_percent": 1.0}) if name.endswith("json") else grid is None
        with open(os.path.join(d, "plain.csv"), "w") as f:                    # no header, no weights
            f.write("1.5, 100\n-1.5, 98\n")
        back, _ = load_points(os.path.join(d, "plain.csv"))
        assert np.allclose(back, [[1500, 1.0, 1.0], [-1500, 0.98, 1.0]])
        with open(os.path.join(d, "bad.json"), "w") as f:
            f.write('{"points": [[0, 100, 1]]}')
        try:
            load_points(os.path.join(d, "bad.json"))
        except ValueError:
            pass
        else:
            raise AssertionError("loaded a file without the format tag")


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
