"""Point selections of the GRAPE goals: points (offset [Hz], B1 scale [1 = nominal], weight, goal).

goal: None (= the gate goal of the GRAPE tab) or the point's own goal {"gate": "UR", "axis", "angle"} /
{"gate": "PP", "m0", "mt"} (e.g. selective pulses: excite inside a band, do nothing outside).

Text syntax (GUI fields): values separated by commas, semicolons or spaces; ``a:b:n`` is n evenly spaced
values from a to b (both included). Example: "-10:10:21, 12, 15". Offsets are given in kHz, B1 in %.
A grid selection is every offset combined with every B1 value (weight 1).
Files: JSON {"format": "psa-grape-points", "version": 1, "points": [[offset_khz, b1_percent, weight, goal], ...],
"grid": {"offset_khz": ..., "b1_percent": ...}}; CSV with the columns offset_khz, b1_percent, weight, goal
(goal text: empty = default, "UR:angle:ax,ay,az" or "PP:m0x,m0y,m0z:mtx,mty,mtz").
"""
from __future__ import annotations

import csv
import json
import re

import numpy as np

FORMAT, VERSION = "psa-grape-points", 1
DEFAULT_GRID = {"offset_khz": 0.1, "b1_percent": 0.5}       # lock-to-grid steps of the visual editor


def parse_values(text, what="values"):
    """Numbers of a field ("a, b", "a:b:n" ranges; separators , ; or spaces), in input order, duplicates
    removed. Raises ValueError with a readable message."""
    out = []
    tokens = [t for t in re.split(r"[,;\s]+", text.strip()) if t]
    if not tokens:
        raise ValueError(f"{what}: enter at least one value")
    for tok in tokens:
        parts = tok.split(":")
        try:
            if len(parts) == 1:
                vals = [float(parts[0])]
            elif len(parts) == 3:
                a, b, n = float(parts[0]), float(parts[1]), int(parts[2])
                if n < 1:
                    raise ValueError
                vals = list(np.linspace(a, b, n)) if n > 1 else [a]
            else:
                raise ValueError
        except ValueError:
            raise ValueError(f"{what}: '{tok}' is neither a number nor a range a:b:n")
        if not all(np.isfinite(vals)):
            raise ValueError(f"{what}: '{tok}' is not finite")
        out += [float(v) for v in vals]
    seen, uniq = set(), []
    for v in out:
        key = round(v, 12)
        if key not in seen:
            seen.add(key)
            uniq.append(v)
    return uniq


def grid_points(offsets_khz, b1_percent):
    """Points [offset_hz, scale, 1] of every offset with every B1 value (offset-major)."""
    if any(b <= 0 for b in b1_percent):
        raise ValueError("B1 values must be positive (in % of the nominal amplitude)")
    return [[o * 1e3, b / 100.0, 1.0] for o in offsets_khz for b in b1_percent]


def is_grid(points):
    """(offsets_khz, b1_percent) if the points are a full grid with all weights 1 and no own goals, else None."""
    if not points or any(abs(p[2] - 1.0) > 1e-12 or (len(p) > 3 and p[3] is not None) for p in points):
        return None
    offs = sorted({round(p[0] / 1e3, 9) for p in points})
    b1 = sorted({round(p[1] * 100, 9) for p in points})
    have = {(round(p[0] / 1e3, 9), round(p[1] * 100, 9)) for p in points}
    if len(have) != len(points) or len(have) != len(offs) * len(b1):
        return None
    return offs, b1


def snap(offset_khz, b1_percent, grid=DEFAULT_GRID):
    """Nearest point of the lock-to-grid lattice."""
    so, sb = float(grid["offset_khz"]), float(grid["b1_percent"])
    o = round(round(offset_khz / so) * so, 9) if so > 0 else offset_khz
    b = round(round(b1_percent / sb) * sb, 9) if sb > 0 else b1_percent
    return o, b


def format_values(vals):
    """Compact text of values (a:b:n for an evenly spaced run of >= 4 values)."""
    vals = [float(v) for v in vals]
    if len(vals) >= 4:
        d = np.diff(vals)
        if abs(d[0]) > 0 and np.allclose(d, d[0], rtol=0, atol=1e-9 * max(1.0, abs(d[0]))):
            return f"{vals[0]:g}:{vals[-1]:g}:{len(vals)}"
    return ", ".join(f"{v:g}" for v in vals)


AXES = {(1, 0, 0): "x", (-1, 0, 0): "−x", (0, 1, 0): "y", (0, -1, 0): "−y", (0, 0, 1): "z", (0, 0, -1): "−z"}


def _vec_text(v):
    u = np.asarray(v, float)
    u = u / np.linalg.norm(u)
    key = tuple(int(round(c)) for c in u)
    if np.allclose(u, key, atol=1e-9) and key in AXES:
        return AXES[key]
    return "(" + ", ".join(f"{c:.3g}" for c in np.asarray(v, float)) + ")"


def describe_goal(goal):
    """Short display text of a goal: "default", "UR 90° x", "PP z→x", ..."""
    if goal is None:
        return "default"
    if goal["gate"] == "UR":
        if abs(float(goal["angle"])) < 1e-12:
            return "UR 0° (nothing)"
        return f"UR {float(goal['angle']):g}° {_vec_text(goal['axis'])}"
    return f"PP {_vec_text(goal['m0'])}→{_vec_text(goal['mt'])}"


def goal_to_text(goal):
    """Compact text of a goal for CSV ("" for the default)."""
    if goal is None:
        return ""
    v = lambda x: ",".join(f"{float(c):.12g}" for c in x)
    if goal["gate"] == "UR":
        return f"UR:{float(goal['angle']):.12g}:{v(goal['axis'])}"
    return f"PP:{v(goal['m0'])}:{v(goal['mt'])}"


def goal_from_text(text):
    """Inverse of goal_to_text (ValueError if invalid)."""
    text = (text or "").strip()
    if not text or text.lower() == "default":
        return None
    parts = text.split(":")
    try:
        vec = lambda t: [float(c) for c in t.split(",")]
        if parts[0].upper() == "UR" and len(parts) == 3:
            goal = {"gate": "UR", "angle": float(parts[1]), "axis": vec(parts[2])}
        elif parts[0].upper() == "PP" and len(parts) == 3:
            goal = {"gate": "PP", "m0": vec(parts[1]), "mt": vec(parts[2])}
        else:
            raise ValueError
    except ValueError:
        raise ValueError(f"invalid goal '{text}' (UR:angle:ax,ay,az or PP:m0x,m0y,m0z:mtx,mty,mtz)")
    from psa.design.grape import check_goal
    check_goal(goal)
    return goal


def save_points(path, points, grid=None):
    """Write a selection (points [offset_hz, scale, weight, goal]) as JSON or, for *.csv, as CSV."""
    rows = [[p[0] / 1e3, p[1] * 100, p[2], p[3] if len(p) > 3 else None] for p in points]
    if str(path).lower().endswith(".csv"):
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["offset_khz", "b1_percent", "weight", "goal"])
            w.writerows([r[:3] + [goal_to_text(r[3])] for r in rows])
    else:
        with open(path, "w") as f:
            json.dump({"format": FORMAT, "version": VERSION, "points": rows, "grid": dict(grid or DEFAULT_GRID)}, f,
                      indent=1)


def load_points(path):
    """(points [offset_hz, scale, weight, goal], grid steps or None) of a JSON or CSV file."""
    if str(path).lower().endswith(".csv"):
        with open(path, newline="") as f:
            rows = [r for r in csv.reader(f) if r and any(c.strip() for c in r)]
        if rows and not _is_number(rows[0][0]):
            rows = rows[1:]
        data, grid = [[float(c) for c in r[:3]] + [1.0] * (3 - len(r[:3])) + [goal_from_text(r[3]) if len(r) > 3 else None]
                      for r in rows], None
    else:
        with open(path) as f:
            d = json.load(f)
        if d.get("format") != FORMAT:
            raise ValueError("not a psa GRAPE point selection")
        data, grid = d["points"], d.get("grid")
    pts = []
    from psa.design.grape import check_goal
    for r in data:
        if len(r) < 2 or float(r[1]) <= 0 or (len(r) > 2 and float(r[2]) <= 0):
            raise ValueError(f"invalid point {r}: B1 and weight must be positive")
        goal = r[3] if len(r) > 3 else None
        if goal is not None:
            check_goal(goal)
        pts.append([float(r[0]) * 1e3, float(r[1]) / 100.0, float(r[2]) if len(r) > 2 else 1.0, goal])
    if not pts:
        raise ValueError("the file has no points")
    return pts, grid


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
