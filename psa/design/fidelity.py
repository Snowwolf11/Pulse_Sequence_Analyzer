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
Freedom (UR / PP / own): the reference may hold only up to a z-rotation after the pulse, theta = phi + 2 pi
offset tau ("phase": tau = 0; "linear": phase + delay, removable by a phase correction / a delay). phi and tau
are either fitted on an area (the map's offset × amplitude grid, as the GRAPE optimizer fits them at its points)
or fixed values (the correction actually applied in the experiment).

A job (``run_job``) computes any of: the map over an offset × amplitude grid, a cut over offsets at a
fixed amplitude, a cut over amplitudes at a fixed offset, and the phase profile of the transverse
magnetization with its linear fit (``phase_profile``). ``fidelity_worker`` runs jobs in a separate
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


def _points_config(ref):
    from psa.design.grape import GrapeConfig
    return GrapeConfig(mode="broadband", gate=ref.get("gate", "none"), target_axis=ref.get("axis", [1, 0, 0]),
                       target_angle_deg=ref.get("angle", 180.0), m0=ref.get("m0", [0, 0, 1]),
                       m_target=ref.get("mt", [1, 0, 0]), points=ref["points"], freedom=ref.get("freedom", "none"),
                       max_delay_us=ref.get("max_delay_us", 0.0))


def resolve(ps: PulseSequence, reference=None):
    """The reference with the pulse-dependent parts filled in: "own" gets its matrix "R"."""
    ref = dict(reference or {"kind": "own"})
    if ref.get("kind", "own") == "own" and "R" not in ref:
        ref["R"] = matrices(ps, [0.0], [1.0])[0, 0]
    return ref


def infidelity_grid(ps: PulseSequence, offsets_hz, scales, reference=None):
    """1 - F at every (scale, offset): array (len(scales), len(offsets_hz)). A UR / PP / own reference may carry
    a z-rotation after the pulse: "phi" [rad] and "tau" [s] (see fit_correction)."""
    ref = dict(reference or {"kind": "own"})
    if ref.get("kind", "own") == "points":   # per-point goals of the GRAPE tab: the nearest point's goal
        from psa.design.grape import goal_infidelity_grid
        return goal_infidelity_grid(ps, _points_config(ref), offsets_hz, scales)
    return infidelity_from_matrices(matrices(ps, offsets_hz, scales), resolve(ps, ref), offsets_hz)


def _target(ref):
    """(Rt, m0, mt) of a UR / own / PP reference (Rt None for PP)."""
    kind = ref.get("kind", "own")
    if kind == "PP":
        return None, _unit(ref.get("m0", (0, 0, 1))), _unit(ref.get("mt", (1, 0, 0)))
    if kind == "UR":
        return rotation_from_axis_angle(_unit(ref.get("axis", (1, 0, 0))), float(ref.get("angle", 180.0))), None, None
    if kind == "own" and "R" in ref:
        return np.asarray(ref["R"], float), None, None
    raise ValueError(f"unknown reference {kind!r}")


def infidelity_from_matrices(R, reference, offsets_hz=None):
    """1 - F of rotation matrices R (..., n_offsets, 3, 3) for a reference dict (UR / PP; "own" with its matrix
    "R"), optionally up to the z-rotation phi + 2 pi offset tau after the pulse ("phi" [rad], "tau" [s])."""
    ref = dict(reference)
    Rt, m0, mt = _target(ref)
    phi, tau = float(ref.get("phi", 0.0)), float(ref.get("tau", 0.0))
    if phi or tau:
        from psa.design.grape import rz
        if tau and offsets_hz is None:
            raise ValueError("a linear phase needs the offsets")
        f = np.zeros(1) if offsets_hz is None else np.asarray(offsets_hz, float)
        Z = rz(phi + 2 * np.pi * f * tau)                                   # (n_offsets, 3, 3)
        Rt, mt = (None, Z @ mt) if Rt is None else (Z @ Rt, None)
    if Rt is None:
        return np.sum((R @ m0 - mt) ** 2, axis=-1) / 4
    return np.sum((R - Rt) ** 2, axis=(-2, -1)) / 12


def fit_correction(ps: PulseSequence, reference, offsets_hz, scales, freedom, max_delay_us=0.0, R=None):
    """(phi [rad], tau [s]) of the best z-rotation after the pulse on the offset × amplitude area (equal weights;
    tau = 0 for "phase"; |tau| <= max_delay_us, 0 = the pulse duration). R: the area's matrices if known."""
    from psa.design.grape import _overlaps, fit_freedom
    ref = resolve(ps, reference)
    Rt, m0, mt = _target(ref)
    if R is None:
        R = matrices(ps, offsets_hz, scales)
    a, b = _overlaps(R, ("PP", None, m0, mt) if Rt is None else ("UR", Rt, None, None))
    f = np.broadcast_to(np.asarray(offsets_hz, float), a.shape)
    tau_max = max_delay_us * 1e-6 if max_delay_us > 0 else ps.duration
    return fit_freedom(a.ravel(), b.ravel(), f.ravel(), np.ones(a.size), freedom, tau_max)


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


def phase_profile(ps: PulseSequence, offsets_hz, scale=1.0, m0=(0, 0, 1), max_delay_us=0.0):
    """Phase of the transverse magnetization R m0 vs. offset (at one amplitude scaling) and its linear fit
    phase(offset) = phi + 360 offset tau, fitted like the "linear phase" freedom (weighted by |m_xy|; |tau| <=
    max_delay_us, 0 = the pulse duration): tau is the delay a first-order phase correction removes.
    Returns {"offsets_hz", "mxy" (|m_xy|), "phase_deg" (unwrapped, next to the line), "fit_deg" (the line),
    "residual_deg" (phase - line, wrapped to ±180), "phi_deg", "tau_us"}."""
    from psa.design.grape import fit_freedom
    f = np.asarray(offsets_hz, float)
    m = matrices(ps, f, [float(scale)])[0] @ _unit(m0)                      # (n_offsets, 3)
    mx, my = m[:, 0], m[:, 1]
    tau_max = max_delay_us * 1e-6 if max_delay_us > 0 else ps.duration
    phi, tau = fit_freedom(mx, my, f, np.ones(len(f)), "linear", tau_max)
    line = phi + 2 * np.pi * f * tau
    resid = np.angle(np.exp(1j * (np.arctan2(my, mx) - line)))
    return {"offsets_hz": f, "mxy": np.hypot(mx, my), "phase_deg": np.degrees(line + resid),
            "fit_deg": np.degrees(line), "residual_deg": np.degrees(resid),
            "phi_deg": float(np.degrees(np.angle(np.exp(1j * phi)))), "tau_us": tau * 1e6}


def run_job(job):
    """job: {"id", "pulse": [amp, phase, dt, umax], "reference", "map": {"offsets_hz", "scales"} | None,
    "cut_offset": {"offsets_hz", "scale"} | None, "cut_amp": {"scales", "offset_hz"} | None,
    "fit_area": {"offsets_hz", "scales"} (where a free correction is fitted; default: the map),
    "phase": {"offsets_hz", "scale", "m0", "max_delay_us"} | None (phase profile, see phase_profile)}.
    A UR / PP / own reference may hold "freedom" ("none" / "phase" / "linear") with either "fit": True (and
    "max_delay_us") or fixed "phi_deg" / "tau_us".
    Returns {"id", "map", "cut_offset", "cut_amp" (infidelity arrays or None), "correction" ({"phi_deg", "tau_us",
    "source"} or None), "elapsed_s"} or {"id", "error"}."""
    t0 = time.time()
    try:
        amp, phase, dt, umax = job["pulse"]
        ps = PulseSequence(np.asarray(amp, float), np.asarray(phase, float), dt, umax)
        ref = dict(job.get("reference") or {"kind": "own"})
        out = {"id": job["id"], "map": None, "cut_offset": None, "cut_amp": None, "correction": None, "phase": None}
        R_map = None
        if ref.get("kind") == "points":
            if ref.get("freedom", "none") != "none":
                from psa.design.grape import fitted_freedom
                phi, tau = fitted_freedom(ps, _points_config(ref))
                out["correction"] = {"phi_deg": np.degrees(phi), "tau_us": tau * 1e6, "source": "fitted at the GRAPE points"}
        elif ref.get("freedom", "none") != "none":
            ref = resolve(ps, ref)
            if ref.get("fit", True):
                area = job.get("fit_area") or job.get("map")
                if not area:
                    raise ValueError("fitting the correction needs an area (the map ranges)")
                if job.get("map") and area is job.get("map"):
                    R_map = matrices(ps, area["offsets_hz"], area["scales"])
                phi, tau = fit_correction(ps, ref, area["offsets_hz"], area["scales"], ref["freedom"],
                                          ref.get("max_delay_us", 0.0), R_map)
                source = "fitted on the map area"
            else:
                phi, tau = np.radians(float(ref.get("phi_deg", 0.0))), float(ref.get("tau_us", 0.0)) * 1e-6
                source = "fixed"
            if ref["freedom"] == "phase":
                tau = 0.0
            ref["phi"], ref["tau"] = phi, tau
            out["correction"] = {"phi_deg": float(np.degrees(phi)), "tau_us": tau * 1e6, "source": source}
        if job.get("map"):
            m = job["map"]
            out["map"] = (infidelity_from_matrices(R_map, ref, m["offsets_hz"]) if R_map is not None
                          else infidelity_grid(ps, m["offsets_hz"], m["scales"], ref))
        if job.get("cut_offset"):
            c = job["cut_offset"]
            out["cut_offset"] = infidelity_grid(ps, c["offsets_hz"], [c["scale"]], ref)[0]
        if job.get("cut_amp"):
            c = job["cut_amp"]
            out["cut_amp"] = infidelity_grid(ps, [c["offset_hz"]], c["scales"], ref)[:, 0]
        if job.get("phase"):
            c = job["phase"]
            out["phase"] = phase_profile(ps, c["offsets_hz"], c.get("scale", 1.0), c.get("m0", (0, 0, 1)),
                                         c.get("max_delay_us", 0.0))
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
