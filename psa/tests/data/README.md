# Test data

- `golden_reference.npz`: 180 reference arrays produced by the original pure-Python backend
  (commit d20cfc5, before the Rust extensions were replaced by `psa/rotations.py`).
  Covers tangents, total rotations, curves (rotation-matrix method), full analysis (both methods),
  single-point qualities, all stability modes and all quality-image combinations.
- `PS_ur.bruker`: UR360 pulse (1000 sub-pulses). `PS_mix.csv`: mixed test sequence incl. zero amplitude
  and a phase ramp. `PS_short.csv`: first 150 sub-pulses of `PS_ur.bruker` (used for stability).
- `quality_dir/`: three BEBOP pulses and the UR360 pulse, used for the quality-image tests.

The reference results intentionally include known bugs (see BUG_TRACKER, e.g. PSA-01, PSA-04);
regenerate this file only when such a bug is fixed on purpose.
