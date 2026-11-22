"""Tests of the editor backend (psa.design): pulse <-> curve conversion, metrics, file IO.

Run with ``python -m pytest psa/tests`` or ``python psa/tests/test_design_core.py``.
"""
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from psa.design.pulse import PulseSequence
from psa.design.conversion import ps_to_curve, ps_to_tangents, tangents_to_ps, curve_to_ps, align_curve
from psa.design import metrics

DT, UMAX = 0.5e-6, 1e4
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1, -1]).astype(complex)


# ------------------------------------------------------------------ helpers ---

def random_pulse(N=400, seed=0, zero_fraction=0.15):
    """Smooth-ish random pulse with zero-amplitude stretches and some phase jumps."""
    rng = np.random.default_rng(seed)
    amp = np.clip(60 + 40 * np.sin(np.linspace(0, 7, N) + rng.uniform(0, 6)) + rng.normal(0, 5, N), 0, 100)
    phase = np.cumsum(rng.normal(0, 4, N)) + np.where(rng.uniform(size=N) < 0.02, 180.0, 0.0).cumsum()
    start = rng.integers(0, N - N // 5)
    amp[start:start + int(zero_fraction * N)] = 0.0
    return PulseSequence(amp, phase, DT, UMAX)


def su2_tangents(ps):
    """Independent reference: tangents U_k^dag sz U_k by exact 2x2 propagation."""
    U = np.eye(2, dtype=complex)
    out = []
    def t(U):
        M = U.conj().T @ SZ @ U
        return np.real([np.trace(M @ s) / 2 for s in (SX, SY, SZ)])
    out.append(t(U))
    for a, p in zip(ps.amplitude, ps.phase):
        w = 2 * np.pi * UMAX * a / 100 * np.array([np.cos(np.radians(p)), np.sin(np.radians(p))])
        n = np.linalg.norm(w)
        th = n * ps.dt
        H = (w[0] * SX + w[1] * SY) / n if n else 0 * SX
        U = (np.cos(th / 2) * np.eye(2) - 1j * np.sin(th / 2) * H) @ U
        out.append(t(U))
    return np.array(out)


def same_pulse(a, b, amp_tol, phase_tol):
    """Amplitudes equal; phases equal modulo 360 wherever the amplitude matters."""
    assert a.N == b.N
    assert np.max(np.abs(a.amplitude - b.amplitude)) < amp_tol, np.max(np.abs(a.amplitude - b.amplitude))
    on = a.amplitude > 1.0
    d = (a.phase[on] - b.phase[on] + 180) % 360 - 180
    assert np.max(np.abs(d)) < phase_tol, np.max(np.abs(d))


# ---------------------------------------------------------- test curves ---

def helix(n=4000, turns=3.0, radius=1.0, pitch=0.6):
    u = np.linspace(0, 2 * np.pi * turns, n)
    return np.c_[radius * np.cos(u), radius * np.sin(u), pitch * u / (2 * np.pi)]


def lissajous_knot(n=6000, nx=3, ny=2, nz=7, px=0.7, py=0.2):
    u = np.linspace(0, 2 * np.pi, n)
    return np.c_[np.cos(nx * u + px), np.cos(ny * u + py), np.cos(nz * u)]


def torus_knot(n=6000, p=2, q=3, R=2.0, r=0.8):
    u = np.linspace(0, 2 * np.pi, n)
    return np.c_[(R + r * np.cos(q * u)) * np.cos(p * u), (R + r * np.cos(q * u)) * np.sin(p * u), r * np.sin(q * u)]


def planar_s_curve(n=4000):
    """Curvature changes sign (inflection points, kappa = 0 crossings)."""
    u = np.linspace(0, 4 * np.pi, n)
    return np.c_[np.zeros(n), 0.8 * np.sin(u), u]


def lines_and_arcs(density=600):
    """G1 curve: straight - arc (xz plane) - straight - arc (yz plane, opposite sense) - straight.
    Contains zero-curvature parts and curvature jumps."""
    pts, pos, tan = [], np.zeros(3), np.array([0.0, 0, 1])
    def line(length):
        nonlocal pos
        s = np.linspace(0, length, int(density * length) + 2)[1:]
        pts.extend(pos + s[:, None] * tan)
        pos = pos + length * tan
    def arc(radius, angle, normal):
        nonlocal pos, tan
        normal = normal / np.linalg.norm(normal)
        c = pos + radius * normal                       # centre
        e1, e2 = -normal, tan
        a = np.linspace(0, angle, int(density * radius * abs(angle)) + 2)[1:]
        sgn = np.sign(angle)
        pts.extend(c + radius * (np.cos(np.abs(a))[:, None] * e1 + np.sin(np.abs(a))[:, None] * e2))
        pos = c + radius * (np.cos(abs(angle)) * e1 + np.sin(abs(angle)) * e2)
        tan = -np.sin(abs(angle)) * e1 + np.cos(abs(angle)) * e2
        tan /= np.linalg.norm(tan)
    pts.append(pos.copy())
    line(1.0)
    arc(0.5, 2.0, np.array([1.0, 0, 0]))
    line(0.8)
    arc(0.4, 1.5, np.cross(tan, np.array([0, 1.0, 0])) if np.linalg.norm(np.cross(tan, [0, 1.0, 0])) > 1e-9 else np.array([0, 1.0, 0]))
    line(0.6)
    arc(0.7, 1.2, np.array([0, 1.0, 0]) - (np.array([0, 1.0, 0]) @ tan) * tan)
    line(0.5)
    return np.array(pts)


CURVES = {"helix": helix, "lissajous": lissajous_knot, "torus_knot": torus_knot,
          "s_curve": planar_s_curve, "lines_and_arcs": lines_and_arcs}


# -------------------------------------------------------------- conversion ---

def test_forward_matches_su2():
    ps = random_pulse(300, seed=1)
    assert np.max(np.abs(ps_to_tangents(ps) - su2_tangents(ps))) < 1e-10
    assert len(ps_to_curve(ps)) == ps.N + 1


def test_tangent_inversion_is_exact():
    """PS -> exact tangents -> PS reproduces the pulse (incl. zero-amplitude parts, phase jumps)."""
    for seed in range(4):
        ps = random_pulse(500, seed=seed)
        back = tangents_to_ps(ps_to_tangents(ps), DT, UMAX, phase_start=ps.phase[0])
        same_pulse(ps, back, amp_tol=1e-7, phase_tol=1e-6)
        assert np.max(np.abs(ps_to_tangents(back) - ps_to_tangents(ps))) < 1e-10


def test_curve_to_ps_reproduces_curves():
    """curve -> PS -> curve for several curve types; the recomputed curve matches the target."""
    for name, make in CURVES.items():
        pts = make()
        ps, info = curve_to_ps(pts, DT, UMAX)
        assert info["amplitude_ok"], (name, info["max_amplitude_percent"])
        assert np.max(ps.amplitude) > 90, name                  # minimal duration uses the full RF
        curve = ps_to_curve(ps)
        L = info["duration"]
        err = np.max(np.linalg.norm(curve - info["target_points"], axis=1)) / L
        assert err < 2e-4, (name, err)
        # the target points are the original curve (aligned and scaled)
        aligned, _ = align_curve(pts)
        assert np.linalg.norm(info["target_points"][-1] - aligned[-1] * info["scale"]) / L < 1e-4, name


def test_curve_to_ps_with_longer_duration():
    pts = torus_knot()
    ps_min, info = curve_to_ps(pts, DT, UMAX)
    ps, info2 = curve_to_ps(pts, DT, UMAX, duration=2.5 * info["duration"])
    assert abs(info2["duration"] - 2.5 * info["duration"]) <= DT
    assert np.max(ps.amplitude) < 0.45 * 100                     # 2.5x longer -> 2.5x less curvature
    assert np.max(np.linalg.norm(ps_to_curve(ps) - info2["target_points"], axis=1)) / info2["duration"] < 2e-4


def test_zero_curvature_parts_give_zero_amplitude():
    pts = lines_and_arcs()
    ps, info = curve_to_ps(pts, DT, UMAX)
    L = info["length"]
    s = np.linspace(0, 1, ps.N)                                  # fraction of the curve
    first_line = s < 0.8 * 1.0 / L                               # first straight piece (length 1.0), away from the junction
    assert np.max(ps.amplitude[first_line]) < 1.0
    assert np.max(ps.amplitude) > 90


def test_inflection_gives_phase_jump():
    ps, _ = curve_to_ps(planar_s_curve(), DT, UMAX)
    on = ps.amplitude > 20
    # planar curve: phase only takes two values (mod 360) that differ by 180 degrees
    ph = np.mod(ps.phase[on] - ps.phase[on][0] + 90, 360) - 90
    assert np.all((np.abs(ph) < 1.0) | (np.abs(ph - 180) < 1.0))
    assert np.any(np.abs(ph - 180) < 1.0)


def test_pulse_curve_pulse_round_trip():
    """PS -> dense curve -> PS -> curve: the pulse (and its curve) are recovered."""
    rng = np.random.default_rng(7)
    N = 300
    amp = 70 + 25 * np.sin(np.linspace(0, 5, N))
    phase = np.cumsum(rng.normal(0, 2, N))
    ps = PulseSequence(amp, phase, DT, UMAX)
    fine = PulseSequence(np.repeat(amp, 10), np.repeat(phase, 10), DT / 10, UMAX)   # same pulse, dense curve
    pts = ps_to_curve(fine)
    back, info = curve_to_ps(pts, DT, UMAX, duration=ps.duration, align=False)
    same_pulse(ps, back, amp_tol=0.5, phase_tol=1.0)
    assert np.max(np.linalg.norm(ps_to_curve(back) - ps_to_curve(ps), axis=1)) / ps.duration < 1e-4


# ----------------------------------------------------------------- metrics ---

def exact_magnus(ps):
    """Independent A1..A3: sequential block-matrix exponentials (no tree reduction)."""
    from scipy.linalg import expm
    X = [np.eye(2, dtype=complex)] + [np.zeros((2, 2), complex)] * 3
    for a, p in zip(ps.amplitude, ps.phase):
        w = 2 * np.pi * UMAX * a / 100
        A0 = -0.5j * w * (np.cos(np.radians(p)) * SX + np.sin(np.radians(p)) * SY)
        Z = np.zeros((2, 2), complex)
        E = expm(np.block([[A0, Z, Z, Z], [-0.5j * SZ, A0, Z, Z], [Z, -0.5j * SZ, A0, Z], [Z, Z, -0.5j * SZ, A0]]) * ps.dt)
        X = list(np.split(E @ np.vstack(X), 4))
    D = [X[0].conj().T @ x for x in X[1:]]
    M = [D[0], D[1] - D[0] @ D[0] / 2, D[2] - (D[0] @ D[1] + D[1] @ D[0]) / 2 + D[0] @ D[0] @ D[0] / 3]
    return [np.real([1j * 2 ** k * np.trace(m @ s) for s in (SX, SY, SZ)]) for k, m in enumerate(M)]


def test_robustness_exact_vs_independent_and_geometric():
    ps = random_pulse(250, seed=3)
    rb = metrics.robustness(ps)
    ref = exact_magnus(ps)
    geo = metrics.robustness_geometric(ps_to_curve(ps))
    for k in range(3):
        A = rb[f"A{k + 1}"]
        assert np.allclose(A, ref[k], rtol=1e-9, atol=1e-12 * ps.duration ** (k + 1)), k
        assert np.linalg.norm(geo[f"A{k + 1}"] - A) < 2e-3 * np.linalg.norm(A) + 1e-6 * ps.duration ** (k + 1), k
    # the trapezoid curve end point approximates A1 = r(T) to O(dt^2)
    assert np.linalg.norm(rb["A1"] - ps_to_curve(ps)[-1]) < 1e-4 * np.linalg.norm(rb["A1"])


def test_robustness_of_closed_4h_pulse():
    """The 4H pulse is first- and second-order robust: A1 = A2 = 0."""
    ps = PulseSequence.load(os.path.join(HERE, "data", "PS_4H_closed.bruker"))
    rb = metrics.robustness(ps)
    assert rb["A1_norm"] < 1e-9 and rb["A2_norm"] < 1e-9
    assert rb["A3_norm"] > 1e-6                                 # third order is not cancelled


def test_gate_of_simple_pulses():
    ps = PulseSequence(np.full(50, 100.0), np.zeros(50), DT, UMAX)          # 25 us at 10 kHz = 90 deg
    g = metrics.gate(ps)
    assert abs(g["angle_deg"] - 90) < 1e-9 and np.allclose(g["axis"], [1, 0, 0])
    ps = PulseSequence(np.full(100, 100.0), np.full(100, 90.0), DT, UMAX)   # 180 deg about y
    g = metrics.gate(ps)
    assert abs(g["angle_deg"] - 180) < 1e-9 and np.allclose(np.abs(g["axis"]), [0, 1, 0], atol=1e-9)
    assert np.allclose(g["matrix"] @ [0, 0, 1], [0, 0, -1], atol=1e-12)


def test_fidelity_vs_offset():
    ps = random_pulse(200, seed=5)
    offs = np.linspace(-5e3, 5e3, 11)
    F = metrics.fidelity_vs_offset(ps, offs)
    assert abs(F[5] - 1) < 1e-12 and np.all(F <= 1 + 1e-12)
    # reference: full SU(2) at 3 kHz
    def U(offset):
        U = np.eye(2, dtype=complex)
        for a, p in zip(ps.amplitude, ps.phase):
            w = 2 * np.pi * np.array([UMAX * a / 100 * np.cos(np.radians(p)), UMAX * a / 100 * np.sin(np.radians(p)), offset])
            n = np.linalg.norm(w)
            if n == 0:
                continue
            H = (w[0] * SX + w[1] * SY + w[2] * SZ) / n
            U = (np.cos(n * ps.dt / 2) * np.eye(2) - 1j * np.sin(n * ps.dt / 2) * H) @ U
        return U
    ref = (abs(np.trace(U(0).conj().T @ U(3e3))) ** 2 + 2) / 6
    assert abs(metrics.fidelity_vs_offset(ps, [3e3])[0] - ref) < 1e-10


# --------------------------------------------------------------------- IO ---

def test_save_and_load():
    ps = random_pulse(120, seed=9)
    ps.dt, ps.umax = 1e-6, 2e4
    with tempfile.TemporaryDirectory() as d:
        ps.save(os.path.join(d, "p.csv"))
        back = PulseSequence.load(os.path.join(d, "p.csv"), dt=1e-6, umax=2e4)
        assert np.allclose(back.as_array(), ps.as_array())
        ps.save(os.path.join(d, "p.bruker"))
        back = PulseSequence.load(os.path.join(d, "p.bruker"))
        assert back.dt == 1e-6 and back.umax == 2e4
        same_pulse(ps, back, amp_tol=1e-9, phase_tol=1e-9)


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
