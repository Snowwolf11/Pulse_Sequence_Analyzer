"""Import of curves from text files: points (x, y, z) or curvature / torsion tables."""
from __future__ import annotations

import math
import re

import numpy as np

from psa.design.conversion import _rotation_from_axis_angle

FORMATS = {
    "points": "curve points: x, y, z per line",
    "s_kappa_tau": "curvature & torsion: s, κ(s), τ(s) per line",
    "kappa_tau": "curvature & torsion: κ, τ per line (equally spaced, step 1)",
}


def load_table(path, min_columns=1):
    """Numeric table from a text / CSV file. Separators: comma, semicolon, tab or spaces. Lines that
    do not start with a number (headers, '#' comments, Bruker '##' lines) are skipped."""
    rows = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            text = line.strip()
            if not text or not re.match(r"^[+\-]?(\d|\.\d)", text):
                continue
            parts = [p for p in re.split(r"[,;\s]+", text) if p]
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    if not rows:
        raise ValueError("no numeric data found")
    width = max(set(len(r) for r in rows), key=lambda w: sum(len(r) == w for r in rows))
    table = np.array([r for r in rows if len(r) == width], dtype=np.float64)
    if width < min_columns:
        raise ValueError(f"expected at least {min_columns} columns, found {width}")
    return table


def curve_from_curvature_torsion(s, kappa, tau):
    """Integrate the Frenet-Serret equations. Curvature and torsion are taken as constant on each
    interval (mean of its end values); the frame (T, N, B) is rotated exactly about the Darboux vector
    tau T + kappa B, positions by the trapezoid rule. Starts at the origin with T = z, N = x, B = y."""
    s = np.asarray(s, dtype=np.float64)
    kappa = np.asarray(kappa, dtype=np.float64)
    tau = np.asarray(tau, dtype=np.float64)
    if not (len(s) == len(kappa) == len(tau)) or len(s) < 2:
        raise ValueError("s, kappa and tau must have the same length (>= 2)")
    if np.any(np.diff(s) <= 0):
        raise ValueError("s must be strictly increasing")
    T, N, B = np.array([0.0, 0, 1]), np.array([1.0, 0, 0]), np.array([0.0, 1, 0])
    pts = [np.zeros(3)]
    for i in range(len(s) - 1):
        ds = s[i + 1] - s[i]
        k = (kappa[i] + kappa[i + 1]) / 2
        t = (tau[i] + tau[i + 1]) / 2
        w = t * T + k * B
        angle = np.linalg.norm(w) * ds
        T_old = T
        if angle > 0:
            Rm = _rotation_from_axis_angle(tuple(w / np.linalg.norm(w)), angle)
            T, N, B = Rm @ T, Rm @ N, Rm @ B
        pts.append(pts[-1] + (T_old + T) / 2 * ds)
    return np.array(pts)


def import_curve(path, fmt="points"):
    """Curve points from a file in one of FORMATS."""
    if fmt == "points":
        table = load_table(path, 3)
        return table[:, :3]
    if fmt == "s_kappa_tau":
        table = load_table(path, 3)
        return curve_from_curvature_torsion(table[:, 0], table[:, 1], table[:, 2])
    if fmt == "kappa_tau":
        table = load_table(path, 2)
        return curve_from_curvature_torsion(np.arange(len(table), dtype=np.float64), table[:, 0], table[:, 1])
    raise ValueError(f"unknown format {fmt}")
