"""GUI-independent backend of the pulse-sequence editor.

Modules
-------
pulse       PulseSequence data model, CSV / Bruker input and output
conversion  pulse sequence <-> error curve (exact discrete inversion)
metrics     gate (rotation axis / angle), robustness vectors A1..A3, fidelity vs. offset
waveform    keyframe model of the visual amplitude / phase editor
families    parametrized curve families (Lissajous, torus knots, ...)
builder     segment-based curve builder (arcs, helices, lines, Bezier)

Conventions (on resonance, noise axis z): sub-pulse k has Rabi frequency
B1_k/100 * umax (Hz) and phase phi_k (deg). The error curve is unit speed, its tangent is
t(t) = U^dag sigma_z U, its length equals the pulse duration (seconds).
"""
