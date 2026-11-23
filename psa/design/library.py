"""File names and bookkeeping of quick-mode pulses and pulse libraries (series of durations).

Single pulse:   pulse_{xxxx}us_{goal}_{range}_{index}.bruker,  e.g. pulse_0550us_UR235-180_off20_B140_2.bruker
Library:        a folder named like a single pulse, with the duration range instead of one duration
                (pulse_{xxxx}-{yyyy}us_{goal}_{range}_{index}), holding pulse_{xxxx}us.bruker per duration and
                library.json (settings and the score of every pulse).

    xxxx    duration in µs, rounded to the nearest µs, 4 digits
    goal    the goal most points have: UR{axis}-{angle} (angle in degrees, rounded) or PP{m0}-{target}; vectors
            as single digits, the smallest integer direction within 0.5° of the vector (else the closest with
            digits <= 9); a negative component is written with "m" (0m10 = (0, -1, 0)), as is a negative angle
    range   off{width of the offset range in kHz}_B1{width of the B1 range in %} over all points (rounded to
            0.1, without a trailing .0, "p" for the point: off7p5; B10 = no B1 range)
    index   1, 2, ...: the smallest one that does not overwrite an existing file or folder
Numbers are rounded half up (0.5 -> 1).
"""
from __future__ import annotations

import itertools
import json
import os
from decimal import ROUND_HALF_UP, Decimal

import numpy as np

LIBRARY_INDEX = "library.json"
FORMAT, VERSION = "psa-pulse-library", 1


def round_half_up(x, places=0):
    """x rounded half up (away from zero for negative x) to ``places`` decimals, as a Decimal."""
    q = Decimal(1).scaleb(-places)
    return Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP)


def number_text(x, places=1):
    """Rounded number for file names: no trailing .0, "p" for the decimal point, "m" for minus: 20 -> "20",
    7.25 -> "7p3", -90 -> "m90"."""
    d = round_half_up(abs(x), places)
    t = format(d.normalize(), "f").replace(".", "p") if d != 0 else "0"
    return ("m" if x < 0 and d != 0 else "") + t


def _directions(max_digit=9):
    vs = np.array([v for v in itertools.product(range(-max_digit, max_digit + 1), repeat=3) if any(v)], float)
    return vs, np.max(np.abs(vs), axis=1)


_DIRS = None


def int_vector(v, tol_deg=0.5, max_digit=9):
    """Integer direction (components -9..9) of a vector: the one with the smallest largest digit within tol_deg,
    else the closest one."""
    global _DIRS
    if _DIRS is None:
        _DIRS = _directions(max_digit)
    vs, size = _DIRS
    u = np.asarray(v, float)
    n = np.linalg.norm(u)
    if n == 0 or not np.all(np.isfinite(u)):
        raise ValueError("the vector must be finite and non-zero")
    cos = (vs @ (u / n)) / np.linalg.norm(vs, axis=1)
    ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
    ok = np.where(ang <= tol_deg)[0]
    if ok.size:
        i = ok[np.lexsort((ang[ok], size[ok]))[0]]               # smallest digits, then the closest
    else:
        i = int(np.argmin(ang))
    return tuple(int(c) for c in vs[i])


def vector_text(v):
    return "".join(("m" if c < 0 else "") + str(abs(c)) for c in int_vector(v))


def goal_tag(goal):
    """UR{axis}-{angle} or PP{m0}-{target} of a goal dict ({"gate": "UR", "axis", "angle"} / {"gate": "PP", ...})."""
    if goal["gate"] == "UR":
        return f"UR{vector_text(goal['axis'])}-{number_text(goal['angle'], 0)}"
    return f"PP{vector_text(goal['m0'])}-{vector_text(goal['mt'])}"


def config_goal(config):
    """The goal most points of a GrapeConfig have (their own goal or the gate goal); None if there is none."""
    from psa.design.grape import goal_key
    gate = ({"gate": "UR", "axis": list(config.target_axis), "angle": float(config.target_angle_deg)}
            if config.gate == "UR" else
            {"gate": "PP", "m0": list(config.m0), "mt": list(config.m_target)} if config.gate == "PP" else None)
    goals = [g if g is not None else gate for g in config.point_goals()]
    counts = {}
    for g in goals:
        if g is not None:
            k = goal_key(g)
            counts[k] = (counts.get(k, (0, g))[0] + 1, counts.get(k, (0, g))[1])
    if not counts:
        return None
    return max(counts.values(), key=lambda c: c[0])[1]       # ties: the first goal met


def range_tag(config):
    pts = np.array(config.ensemble(), float)
    off = (np.max(pts[:, 0]) - np.min(pts[:, 0])) / 1e3
    b1 = (np.max(pts[:, 1]) - np.min(pts[:, 1])) * 100
    return f"off{number_text(off)}_B1{number_text(b1)}"


def duration_text(duration_s):
    return f"{int(round_half_up(duration_s * 1e6)):04d}us"


def stem(config, duration_text_):
    goal = config_goal(config)
    return f"pulse_{duration_text_}_{goal_tag(goal) if goal is not None else 'nogoal'}_{range_tag(config)}"


def free_path(folder, base, ext=""):
    """folder/base_{i}{ext} with the smallest i >= 1 that does not exist yet."""
    i = 1
    while os.path.exists(os.path.join(folder, f"{base}_{i}{ext}")):
        i += 1
    return os.path.join(folder, f"{base}_{i}{ext}")


def single_path(folder, config, duration_s):
    """Path of a single quick-mode pulse (see the module doc)."""
    return free_path(folder, stem(config, duration_text(duration_s)), ".bruker")


def library_folder(folder, config, durations_s):
    """Path of a new library folder for a series of durations (not created)."""
    d = f"{duration_text(min(durations_s))[:-2]}-{duration_text(max(durations_s))}"
    return free_path(folder, stem(config, d))


def library_file(duration_s):
    return f"pulse_{duration_text(duration_s)}.bruker"


def series_durations(lo_us, hi_us, n):
    """n durations from lo to hi µs (both included), each rounded to the nearest µs; duplicates removed."""
    vals = np.linspace(lo_us, hi_us, int(n)) if int(n) > 1 else np.array([lo_us])
    out = []
    for v in vals:
        r = int(round_half_up(v))
        if r not in out:
            out.append(r)
    return [r * 1e-6 for r in out]


def series_budgets(durations_s, total_s):
    """Time budget of every pulse of a series: the total split in proportion to the durations."""
    d = np.asarray(durations_s, float)
    return list(total_s * d / d.sum())


def write_index(folder, settings, pulses):
    """library.json: the settings of the series and one entry per finished pulse."""
    with open(os.path.join(folder, LIBRARY_INDEX), "w") as f:
        json.dump({"format": FORMAT, "version": VERSION, "settings": settings, "pulses": pulses}, f, indent=1)


def read_index(folder):
    with open(os.path.join(folder, LIBRARY_INDEX)) as f:
        d = json.load(f)
    if d.get("format") != FORMAT:
        raise ValueError("not a psa pulse library")
    return d
