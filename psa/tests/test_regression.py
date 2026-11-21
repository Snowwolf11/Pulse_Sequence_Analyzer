"""Regression test: the NumPy backend must reproduce the reference results of the original
Python implementation (captured before the Rust code was removed, see data/README.md).

Run with ``python -m pytest psa/tests`` or ``python psa/tests/test_regression.py``.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

import psa.app_backend as ab
from psa.createCurve import createCurve
from psa.getPulseSequence import getPulseSequence
from psa.rotations import tangent_vectors, inverse_total_rotation
from psa.calculate_pulse_sequence_quality import calculate_pulse_sequence_quality
from psa.calculate_curve_stability import calculate_curve_stability
from psa.calculate_pulse_sequence_quality_images import calculate_pulse_sequence_quality_images

DATA = os.path.join(HERE, "data")
GOLDEN = np.load(os.path.join(DATA, "golden_reference.npz"))
PS_UR = getPulseSequence(os.path.join(DATA, "PS_ur.bruker"))
PS_MIX = np.loadtxt(os.path.join(DATA, "PS_mix.csv"), delimiter=",")
PS_SHORT = np.loadtxt(os.path.join(DATA, "PS_short.csv"), delimiter=",")
QDIR = os.path.join(DATA, "quality_dir")
VECS = {"z": np.array([0, 0, 1.0]), "x": np.array([1.0, 0, 0]), "d": np.array([1.0, 2.0, -0.5])}
RTOL, ATOL = 1e-9, 1e-12


def close(key, value, rtol=RTOL, atol=ATOL):
    ref = GOLDEN[key]
    value = np.asarray(value)
    assert value.shape == ref.shape, f"{key}: shape {value.shape} != {ref.shape}"
    if ref.dtype.kind in "fc":
        finite = np.isfinite(ref)
        scale = np.max(np.abs(ref[finite])) if finite.any() else 1.0
        # nan/inf entries (e.g. helix-method torsion, see BUG_TRACKER PSA-08) must match exactly
        assert np.allclose(value, ref, rtol=rtol, atol=atol + rtol * scale, equal_nan=True), \
            f"{key}: max abs diff {np.nanmax(np.abs(value - ref)):.3e} (scale {scale:.3e})"
    else:
        assert np.array_equal(value, ref), f"{key}: {value} != {ref}"


def test_tangents_and_total_rotation():
    for pname, PS in (("ur", PS_UR), ("mix", PS_MIX)):
        for off in (0.0, 3000.0):
            for inpo in (1, 4):
                PSi = np.kron(PS, np.ones((inpo, 1))) / [inpo, 1]
                for vn, v in VECS.items():
                    close(f"vec_{pname}_{off:.0f}_{inpo}_{vn}", tangent_vectors(PSi, 0.5e-6, 0.7, 1e4, off, inpo, v / np.linalg.norm(v)))
                close(f"rtot_{pname}_{off:.0f}_{inpo}", inverse_total_rotation(PSi, 0.5e-6, 1e4, off, inpo))


def test_create_curve():
    for pname, PS in (("ur", PS_UR), ("mix", PS_MIX)):
        for off in (0.0, 3000.0):
            for inpo in (1, 5):
                for vn, v in VECS.items():
                    CM, VM, PSm = createCurve(PS.copy(), T=0.5e-6, l=1.3, maximumAmplitude=1e4, offset=off, inpoFact=inpo, calculationMethod=1, initialVector=v)
                    close(f"cm_{pname}_{off:.0f}_{inpo}_{vn}", CM)
                    close(f"vm_{pname}_{off:.0f}_{inpo}_{vn}", VM)
                    close(f"psm_{pname}_{off:.0f}_{inpo}_{vn}", PSm)


def test_analysis():
    for method in (1, 2):
        r = ab.analyze_pulse_sequence(ab.AnalyzerParams(os.path.join(DATA, "PS_ur.bruker"), 10, 1, 0.5, 5, 0, 2, 95, method, np.array([0, 0, 1])))
        for k, v in vars(r).items():
            v = np.asarray(v, dtype=float)
            if method == 2 and k in ("curvature", "torsion", "integrated_curvature", "integrated_torsion",
                                     "integrated_absolut_torsion", "avg_curvature", "avg_torsion"):
                # Helix method duplicates junction points (BUG_TRACKER PSA-08): zero-length segments
                # give inf or ~1e32 depending on NumPy round-off. Compare only the regular points.
                if v.ndim:
                    ref = GOLDEN[f"an2_{k}"]
                    ok = np.isfinite(ref) & np.isfinite(v) & (np.abs(ref) < 1e12) & (np.abs(v) < 1e12)
                    assert np.allclose(v[ok], ref[ok], rtol=1e-7, atol=1e-7 * np.max(np.abs(ref[ok]))), k
                continue
            close(f"an{method}_{k}", v, rtol=1e-7)


def test_single_point_quality():
    for ct in (1, 2, 3):
        for off in (0.0, 2500.0):
            q = calculate_pulse_sequence_quality(PS_MIX, 0.5e-6, 1.0, 9e3, off, 1, np.array([0, 0, 1.0]), np.array([0.6, 0, 0.8]), ct,
                                                 target_Rotation=inverse_total_rotation(PS_MIX, 0.5e-6, 1e4, 0, 1))
            close(f"q_{ct}_{off:.0f}", q)


def test_stability():
    for m in (1, 2, 3):
        X, Y, Z, q, ax, ang = calculate_curve_stability(PS_SHORT, 0.5e-6, 1.0, 1e4, 20, 20, m, np.array([0, 0, 1.0]))
        close(f"st{m}_X", X); close(f"st{m}_Y", Y); close(f"st{m}_Z", Z)
        close(f"st{m}_q", np.atleast_1d(q))
        assert str(ax) == str(GOLDEN[f"st{m}_axis"]) and float(ang) == float(GOLDEN[f"st{m}_angle"])


def test_quality_images():
    for res in (1.0, 2.5):
        for cv in (1, 2):
            for ct in (1, 2, 3):
                PI = calculate_pulse_sequence_quality_images(QDIR, 20000, 0.5e-6, 1e4, np.array([0, 0, 1.0]), ct, cv, res)
                # log10 of infidelities near 1e-15 amplifies round-off; compare in absolute log units
                close(f"qi_{res}_{cv}_{ct}", PI, rtol=0, atol=1e-6)


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
