# Test data

- `golden_reference.npz`: 180 reference arrays produced by the original pure-Python backend
  (commit d20cfc5, before the Rust extensions were replaced by `psa/rotations.py`).
  Covers tangents, total rotations, curves (rotation-matrix method), full analysis (both methods),
  single-point qualities, all stability modes and all quality-image combinations.
- `PS_ur.bruker`: UR360 pulse (1000 sub-pulses). `PS_mix.csv`: mixed test sequence incl. zero amplitude
  and a phase ramp. `PS_short.csv`: first 150 sub-pulses of `PS_ur.bruker` (used for stability).
- `quality_dir/`: three BEBOP pulses and the UR360 pulse, used for the quality-image tests.

Updated 2026-11-21 for the PSA-01 fix (PP fidelity): entries q_1_*, st1_Z, st1_q were recomputed and
checked against exact SU(2) propagation; all other entries are unchanged.
Updated 2026-11-21 for the PSA-04 fix (trapezoid rule, true tangents): cm_*, an1_* (curve-derived),
an2_VM/totalRot/phi, q_3_*, qi_*_3 recomputed and checked against exact SU(2) propagation.
`PS_4H_closed.bruker`: 4H pulse with an exactly closed error curve (4Hpulses/optima1/pulse0400).

The reference results intentionally include known bugs (see BUG_TRACKER, e.g. PSA-04);
regenerate this file only when such a bug is fixed on purpose.
