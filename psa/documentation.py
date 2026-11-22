DOCUMENTATION_STR = """Documentation

Overview
This app computes and visualizes error curves for NMR pulse sequences using the space-curve quantum control formalism. It can also evaluate fidelities and related metrics.

Main Window
- Load pulse sequence:
  • Enter the full file path (e.g., /Users/leon/Documents/MATLAB/Pulssequenzen/BIBOP_sorted_20kHz_noB1_rf10kHz/pulse0015.bruker), or click **Browse**.
- Pulse parameters:
  • Set sub-pulse duration and other parameters in the fields on the left.
- Controls:
  • Press **Refresh** (or hit **Enter**) to recompute the Error Curve, Curvature, Torsion, and Path (trajectory of the vector tip).
  • Adjust **Offset** (kHz) and **Amplitude** (%) via sliders or their input fields.
- Displays:
  • **Curve Data**: numerical summaries for the current sequence and its error curve.
  • **Info & Errors**: status messages and diagnostics.
- Shortcuts:
  • **Edit** opens the Pulse Sequence Editor.
  • **Menu** opens global settings.

Edit Window
- Tkinter / Dash interface: define a pulse sequence as a list of rectangular pulses and phase ramps, using the format shown in the editor when it opens. Choose a target directory; the defined sequence will be saved there. Use **Plot** to visualize the sequence in the Curve Window. Calculations use the pulse parameters from the Main Window.
- PySide6 interface (python main.py --gui pyside6): **Edit** opens the Pulse Sequence Editor with visual editing, curve families, a curve builder, curve import and a solver. It has its own **Help** button; the full description is in the section "Pulse Sequence Editor" below.

Menu / Settings
- Calculation method: **Rotation Matrix** or **Helix**.
- Initial vector: the **error direction** (noise axis) used to build the error curve.
- Initial magnetization: the Bloch vector before the pulse, used for the **PP/SS fidelity** (final magnetization compared with its on-resonance value) and for the **Angle** stability mode (flip angle of the final magnetization from +z). Default: z.

Error-Curve Window
Plots the Error curve. Updates when you press **Refresh**, hit **Enter**, or change Offset/Amplitude.

Trajectory Window
Plots the trajectory of the vector tip (the first derivative of the curve). Updates on **Refresh**, **Enter**, or Offset/Amplitude changes.

Curvature Window
Plots curvature along the curve. Updates on **Refresh**, **Enter**, or Offset/Amplitude changes.

Torsion Window
Plots torsion along the curve. Updates on **Refresh**, **Enter**, or Offset/Amplitude changes.
Note: torsion is only accurate when `inpoFact = 1`.

Stability Window
- Choose pulse type: **State-to-State (SS)** or **Universal Rotation (UR)**.
- Set ranges:
  • **Offset range** in kHz (0 kHz is center; e.g., 20 kHz → −10 kHz … +10 kHz).
  • **Amplitude range** in % (100% is center; e.g., 20% → 90% … 110%).
- Click **Calculate** to plot fidelity (or rotation angle if **Angle** is selected) vs. offset and amplitude.
- PP/SS uses the initial magnetization from the settings; the target is its on-resonance final state.

PS Window
Shows the pulse sequence as plain text.

Options Window

1) Export Data
- **Export Data**: writes analysis for the sequence specified to the right of the button into a newly created directory:
  • CSV: curve, first derivative, curvature, torsion
  • TXT: additional metadata and settings
- **Export Dir Data**: performs the same export for all pulse files in the sequence’s directory, creating one subdirectory per sequence.

2) Create Quality Images
Generates a heatmap image of infidelity vs. offset for sequences of increasing length.
- **Input directory**: folder containing the pulse sequences.
- **Output directory**: where to save the generated image.
- **Offset Range** (kHz): evaluation window (0 kHz is center; e.g., 20 kHz → −10 kHz … +10 kHz).
- **Class**: UR or SS.
- **Resolution**: pixels per kHz per pulse.

Tips
- Press **Enter** in any numeric field to recompute quickly.

"""


EDITOR_DOCUMENTATION_STR = """Pulse Sequence Editor

Overview
The editor creates and optimizes pulse sequences in several ways (tabs on the left) and evaluates them live
(panel on the right). All tabs produce the same thing: N sub-pulses of equal duration with an amplitude
(% of the maximum amplitude) and a phase (degrees). Everything is on resonance; the detuning noise acts
along z. The error curve is unit speed: its length equals the pulse duration, its curvature is the Rabi
frequency (rad/s) and its torsion is the phase change rate. Tabs and panel scroll when the window is small.

Shared panel (right)
- Metadata: sub-pulse duration and maximum amplitude, used by all tabs.
- Resulting pulse: number of sub-pulses and duration, peak amplitude (red if above 100 %), the gate
  (rotation angle and axis of the whole pulse) and its robustness:
  • 1st order |A1|/T: end point of the error curve (0 = closed curve, first-order robust)
  • 2nd order |A2|/T²: vector area of the curve (0 = second-order robust)
  • 3rd order |A3|/T³: third-order term
  • B1 scaling (1st): first-order error for amplitude scaling, |B1| / ∫Ω dt
  Values below 1e-6 are shown in green. All values are exact for the sub-pulse sequence (Magnus terms).
- 3D view of the error curve (start green, end red) and the gate fidelity vs. offset relative to the
  on-resonance gate ("live" can be switched off for long pulses).
- Output: file path and format (CSV or Bruker shape file; the extension decides). "Save + plot in main
  window" saves the pulse and loads it in the analyzer.

Text tab
The pulse syntax R(B1 %, phase °, duration µs); PR(B1 %, start phase °, end phase °, duration µs);
durations must be multiples of the sub-pulse duration. The pulse is parsed while you type.

Waveform tab
Amplitude and phase plots over time with keyframes.
- Drag a keyframe to move it (it snaps to sub-pulse starts), double-click to add one, right-click to delete.
- Click a keyframe to select it: type its exact value, choose "linear to next" or "hold until next".
- Total duration in ms; the phase axis is never wrapped.
- "From Text tab" / "Open pulse file…" import a pulse as keyframes (within the tolerances; 0 = exact).
- Axis ranges: set time (both plots), amplitude and phase ranges, or "Auto". Scrolling over an axis zooms
  only that axis.
- Designs (keyframes) can be saved and opened as JSON.

Curve families tab
Choose a family and change its parameters (slider or number). Curve families (helix, Lissajous knot,
torus knot, Clelia curve, custom formula) are converted exactly into a pulse; the duration is the shortest
possible (full amplitude) or fixed. Pulse families (2H, 4H) are generated directly from their closed-form
construction. "calculate pulse sequence live" can be switched off for slow families: then only the curve
is drawn until you click "Calculate now". The custom family accepts formulas in u with parameters a–d.

Curve builder tab
Build the error curve from segments; each starts where the previous one ended, in the same direction.
Lengths are in µs (= time). R_min = 1/(2π ν_max) is the tightest radius at full amplitude.
- line (length): a pause (amplitude 0).
- arc (radius R ≥ R_min, angle): constant amplitude R_min/R.
- helix (radius ρ, pitch h per turn, angle): a phase ramp; curvature ρ/(ρ²+c²), torsion c/(ρ²+c²), c = h/2π.
- Bézier (control points P1–P3 in the local frame x', y', t): P1 lies on the tangent; drag P1–P3 in the
  projections.
- twist: phase jump at the start of the segment = rotation of the segment about the tangent.
The 3D view and the xy / xz / yz projections show the curve; the selected segment is orange.

Import curve tab
Load curve points (x, y, z per line) or a curvature/torsion table (s, κ, τ or κ, τ at equal spacing)
from a text / CSV file. Units do not matter: the curve is aligned (start at the origin, initial tangent z)
and scaled to the pulse duration.

Solver (Waveform, Curve families: below the plots; Curve builder: left)
Adjusts free variables so that the chosen goals are met:
- cancel detuning error of 1st (closed curve), 2nd (areas) and/or 3rd order,
- B1 scaling robust (1st order),
- optional target gate: rotation angle about an axis.
Free variables: keyframe values (Waveform; their times stay fixed), continuous family parameters
(Curve families; curve families keep the current duration), or the segment parameters ticked "free"
(Curve builder). The solver runs in the background; "Stop" returns the best result found so far.
The result line shows all goals before and after ("solved" when all are below the thresholds of the
solver settings, by default 1e-8 and a gate error of 1e-5°). Not every combination of goals has a
solution for given free variables.
The solver stops (keeping the best result) when
- all goals are reached ("stop as soon as all goals are reached", on by default),
- the evaluation limit is reached (all evaluations counted, as shown in the solver box: "N / limit"),
- the time limit is reached (if set),
- it converges (step, improvement or gradient below xtol / ftol / gtol; 0 switches a criterion off),
- or you click Stop.
Solver settings ("Settings…" in the solver box, or Menu → Settings → Solver in the analyzer): these
stopping rules, the finite-difference step, the "solved" thresholds and a weight per goal (relative
importance in the least-squares residual).

Optimize (GRAPE) tab
Optimizes every sub-pulse individually (Cartesian controls u_x, u_y per sub-pulse) with exact gradients,
like classical GRAPE. The duration and the number of sub-pulses stay those of the start pulse.
- Start pulse: the current pulse of another tab, a rectangular pulse (duration, amplitude, phase) or a
  random one (duration, seed); "Load start pulse". "Result → start" continues from the result.
- Robustness, one of
  • local (SCQC): cancel the detuning error terms A1/A2/A3 (orders ticked) and optionally the first-order
    B1 error, expanded around every offset of the list (kHz, comma-separated, e.g. "-3, 0, 3"; the gate
    goal is then required at every listed offset as well),
  • broadband (classical GRAPE): the gate goal on a grid of offsets × B1 scalings (mean infidelity).
- Gate goal: none, UR (rotation angle about an axis; the whole rotation) or PP (point to point: initial
  magnetization m0 to a target). For PP the robustness terms only count the error that moves m0.
- The amplitude limit (100 %) is always kept exactly (built into the parametrization).
- Controls: "every sub-pulse" (classical) or "smooth spline" with M control points per channel (cubic
  B-splines: smooth by construction, few variables, often much faster; needs M ≤ number of sub-pulses).
- Smoothness weight (0 = off; try 1e-4): prefers smooth u_x, u_y (no spikes). It is soft: the optimizer
  first solves with it, then lowers it step by step (×0.1) until the goals are reached exactly, so the
  result stays close to the smooth solution.
- Max. amplitude rate (%/µs) and max. phase rate (kHz = torsion of the error curve / 2π): soft limits that
  cost nothing below the limit (kept to about 1 %). The phase rate counts in proportion to the amplitude
  (a phase jump where the amplitude is ~0 is harmless). Too strict limits make the goals unreachable.
- Regularization "stay close to start": adds λ·|u − u_start|² (also lowered step by step).
- Starts / spread (exploration): several optimizations from perturbed starts run in parallel processes
  (one residual line each); the result is the best one: goals reached first, then the smoothest (if a
  smoothness weight is set) or the lowest residual. Spread 0.5 is a large perturbation, 0.05 a small one.
- Start noise: a tiny random change of the first start (seeded). An exactly symmetric start, e.g. a
  constant phase, can otherwise keep the optimizer in a symmetric family of pulses that cannot reach the goal.
The optimizer (Levenberg–Marquardt with exact Jacobians) runs in the background (Optimize / Stop); the plots
show the start (dashed), the best pulse so far and the goal residual vs. evaluations. Stopping rules (the
evaluation limit counts per start), "solved" thresholds and goal weights are the shared solver settings
("Settings…"; the finite-difference step is not used: the gradients are exact).
If all starts converge to about the same residual without reaching the goals, the pulse is most likely too
short for them (or the limits too strict): try a longer start pulse.
"""

# the main window's documentation includes the editor description
DOCUMENTATION_STR = DOCUMENTATION_STR.rstrip() + "\n\n\n" + EDITOR_DOCUMENTATION_STR
