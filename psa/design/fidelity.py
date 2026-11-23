"""Fidelity of a pulse over offset × amplitude (B1) scaling, for the editor's fidelity analysis window.

Reference ("what should the pulse do"):
    own   the pulse's own gate at offset 0 and 100 % amplitude (as the editor panel's plot)
    UR    a rotation by an angle about an axis (the whole rotation)
    PP    point to point: initial magnetization m0 to a target
Infidelity 1 - F, computed from the residual (no cancellation for small values):
    UR / own   1 - (3 + tr(Rt^T R)) / 6 = |R - Rt|_F^2 / 12   (= the gate fidelity (|tr U0^† U|^2 + 2) / 6)
    points     the goals of the GRAPE tab at its points (per-point goals, freedom): every location is judged by
               the goal of its nearest point (NaN where that point has no goal); with freedom (free phase /
               linear phase) the targets are rotated by the z-rotation fitted at the points
    PP         1 - (1 + (R m0) . mt) / 2 = |R m0 - mt|^2 / 4

A job (``run_job``) computes any of: the map over an offset × amplitude grid, a cut over offsets at a
fixed amplitude, and a cut over amplitudes at a fixed offset. ``fidelity_worker`` runs jobs in a separate
process (the editor's window does not block, and it never touches running optimizations).
"""
from __future__ import annotations

import time
import traceback

import numpy as np

from psa.design.optimize import rotation_from_axis_angle
from psa.design.pulse import PulseSequence


def _unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("the reference vectors must not be zero")
    return v / n


def matrices(ps, offsets_hz, scales):
    """Rotation matrices of the pulse at every (scale, offset): shape (n_scales, n_offsets, 3, 3)."""
    from psa.rotations import propagate_batch, quaternion_to_matrix
    sc, of = np.meshgrid(np.asarray(scales, float), np.asarray(offsets_hz, float), indexing="ij")
    q = propagate_batch(ps.as_array(), ps.dt, ps.umax * sc, of)
    return quaternion_to_matrix(q.reshape(-1, 4)).reshape(sc.shape + (3, 3))


def infidelity_grid(ps: PulseSequence, offsets_hz, scales, reference=None):
    """1 - F at every (scale, offset): array (len(scales), len(offsets_hz))."""
    ref = dict(reference or {"kind": "own"})
    kind = ref.get("kind", "own")
    if kind == "points":                     # per-point goals of the GRAPE tab: the nearest point's goal
        from psa.design.grape import GrapeConfig, goal_infidelity_grid
        cfg = GrapeConfig(mode="broadband", gate=ref.get("gate", "none"), target_axis=ref.get("axis", [1, 0, 0]),
                          target_angle_deg=ref.get("angle", 180.0), m0=ref.get("m0", [0, 0, 1]),
                          m_target=ref.get("mt", [1, 0, 0]), points=ref["points"], freedom=ref.get("freedom", "none"),
                          max_delay_us=ref.get("max_delay_us", 0.0))
        return goal_infidelity_grid(ps, cfg, offsets_hz, scales)
    R = matrices(ps, offsets_hz, scales)
    if kind == "own":
        ref = {"kind": "own", "R": matrices(ps, [0.0], [1.0])[0, 0]}
    return infidelity_from_matrices(R, ref)


def infidelity_from_matrices(R, reference):
    """1 - F of rotation matrices R (..., 3, 3) for a reference dict (UR / PP; "own" with its matrix "R")."""
    ref = dict(reference)
    kind = ref.get("kind", "own")
    if kind == "PP":
        m = R @ _unit(ref.get("m0", (0, 0, 1)))
        return np.sum((m - _unit(ref.get("mt", (1, 0, 0)))) ** 2, axis=-1) / 4
    if kind == "UR":
        Rt = rotation_from_axis_angle(_unit(ref.get("axis", (1, 0, 0))), float(ref.get("angle", 180.0)))
    elif kind == "own" and "R" in ref:
        Rt = np.asarray(ref["R"], float)
    else:
        raise ValueError(f"unknown reference {kind!r}")
    return np.sum((R - Rt) ** 2, axis=(-2, -1)) / 12


def summary(inf, thresholds=(1e-3, 1e-4, 1e-5)):
    """Worst / mean infidelity and the fraction of the values below each threshold."""
    inf = np.asarray(inf, float)
    inf = inf[np.isfinite(inf)]                 # NaN: no goal there
    if inf.size == 0:
        return {"worst": float("nan"), "mean": float("nan"), **{f"below {t:g}": float("nan") for t in thresholds}}
    out = {"worst": float(np.max(inf)), "mean": float(np.mean(inf))}
    for t in thresholds:
        out[f"below {t:g}"] = float(np.mean(inf <= t))
    return out


def run_job(job):
    """job: {"id", "pulse": [amp, phase, dt, umax], "reference", "map": {"offsets_hz", "scales"} | None,
    "cut_offset": {"offsets_hz", "scale"} | None, "cut_amp": {"scales", "offset_hz"} | None}.
    Returns {"id", "map", "cut_offset", "cut_amp" (infidelity arrays or None), "elapsed_s"} or {"id", "error"}."""
    t0 = time.time()
    try:
        amp, phase, dt, umax = job["pulse"]
        ps = PulseSequence(np.asarray(amp, float), np.asarray(phase, float), dt, umax)
        ref = job.get("reference")
        out = {"id": job["id"], "map": None, "cut_offset": None, "cut_amp": None}
        if job.get("map"):
            m = job["map"]
            out["map"] = infidelity_grid(ps, m["offsets_hz"], m["scales"], ref)
        if job.get("cut_offset"):
            c = job["cut_offset"]
            out["cut_offset"] = infidelity_grid(ps, c["offsets_hz"], [c["scale"]], ref)[0]
        if job.get("cut_amp"):
            c = job["cut_amp"]
            out["cut_amp"] = infidelity_grid(ps, [c["offset_hz"]], c["scales"], ref)[:, 0]
        out["elapsed_s"] = time.time() - t0
        return out
    except Exception as e:
        return {"id": job.get("id"), "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc(limit=3)}


def fidelity_worker(jobs, results):
    """Worker process: runs jobs from ``jobs`` until None; puts the results into ``results``."""
    while True:
        job = jobs.get()
        if job is None:
            return
        results.put(run_job(job))
