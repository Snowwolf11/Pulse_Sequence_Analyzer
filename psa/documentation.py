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
  • **Curve Data** (PySide6): the pulse (sub-pulses, duration, max. amplitude, average power in % of a
    rectangular pulse and as RMS amplitude), its gate on resonance (angle, axis), the robustness |A1|/T,
    |A2|/T², |A3|/T³ and the first-order B1 scaling (as the editor's panel), and the shown curve at the set
    offset / amplitude: the rotation there, length, end point r(T), vector area ½∫r×dr (exact for the
    polygon), integrated curvature / torsion and how far the initial vector turns.
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
Note: torsion is only accurate when `inpoFact = 1` (the default).

Stability Window (PySide6)
The fidelity analysis of the editor (see "Fidelity analysis window" below): fidelity over offset ×
amplitude as 3D surface or heat map with ranges and resolution, cuts at a fixed amplitude / offset,
reference the pulse's own gate (default) / UR / PP, optional free or linear phase; computed on Refresh in
a separate process. Only the main window's pulse (nominal amplitude); no GRAPE goal here.

Stability Window (Tkinter / Dash)
- Choose pulse type: **State-to-State (SS)** or **Universal Rotation (UR)**.
- Set ranges:
  • **Offset range** in kHz (0 kHz is center; e.g., 20 kHz → −10 kHz … +10 kHz).
  • **Amplitude range** in % (100% is center; e.g., 20% → 90% … 110%).
- Click **Calculate** to plot fidelity (or rotation angle if **Angle** is selected) vs. offset and amplitude.
- PP/SS uses the initial magnetization from the settings; the target is its on-resonance final state.

Animation Window (PySide6, needs OpenGL: python -m pip install PyOpenGL)
The spin dynamics of the main window's pulse: the three Bloch-sphere panels in the rotating frame (rotating
at the carrier frequency: an on-resonance spin feels only B1), the error curve in the toggling frame (the
frame that follows the nominal rotation, where only the error accumulates). Four panels side by side (each can
be hidden; the three Bloch-sphere views rotate together, drag to rotate, wheel to zoom):
  • Ensemble: every spin of the ensemble (offsets in kHz and B1 values in %, every offset with every
    value) on the Bloch sphere, colour = offset (blue → red, colour bar), size = B1, fading trails; the
    selected spin is ringed; optional target direction (yellow ring).
  • Effective field: the selected spin's magnetization (cyan), the effective field of the running
    sub-pulse (orange: B1 and offset) and the precession circle around it.
  • Rotation: the selected spin's rotation R(t) as the rotated frame x / y / z (red / green / blue, the final
    frame faint), as DROPS droplets of the propagator U(t) (top rank 0 = identity part, bottom rank 1 =
    rotation-axis part; colour = phase: red 0, yellow 90°, green 180°, blue 270°), or as the Q-bead of the
    magnetization (sphere coloured red along m, green opposite).
  • Error curve: the selected spin's error curve growing with time (green → red) with its tangent.
Below: the pulse with a draggable time cursor, Play / Pause, time slider, playback duration, loop, trail
length, start m0 (the magnetization of every spin at t = 0, e.g. z or x), the selected spin (offset, B1)
and Export GIF (the visible panels, one pass of the pulse).
The motion is exact at every instant (within a sub-pulse the rotation axis is constant).

PS Window
Shows the pulse sequence as plain text (PySide6: sub-tab "Text", one line per sub-pulse with start time,
amplitude and phase; sub-tab "Plots": amplitude, phase (−180 … 180°) and the x / y amplitude components
A cos φ, A sin φ vs. time).

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
- ⊕ (next to the fidelity plot) opens the fidelity analysis window: fidelity vs. offset and amplitude
  (3D surface or heat map; offset and amplitude range and number of points settable), and below it
  fidelity vs. offset at a fixed amplitude and vs. amplitude at a fixed offset, and at the bottom the phase
  of the transverse magnetization (of m0 = z, or the PP reference's m0) vs. offset with its linear fit
  phase + 360° · offset · delay (fitted as the "linear phase" freedom, weighted by the transverse
  magnitude; the delay is the first-order phase correction) and the deviation from that line (rms / max
  where |m_xy| ≥ 0.5). Reference: automatic
  (default: the gate goal of the GRAPE tab while it optimizes or is shown, else the pulse's own gate at
  0 kHz / 100 % as in the panel), or fixed: the own gate, the GRAPE gate goal, or a UR rotation (angle,
  axis) / PP transfer (m0 → target) set in the window (prefilled from the GRAPE goal). The line below the
  choice and the status after Refresh name the reference used. Freedom (own / UR / PP reference): judge
  the pulse only up to a z-rotation after it, a free phase or a linear phase (phase + delay, removable by a
  first-order phase correction or a delay), either fitted on the map's offset × amplitude range (the best
  correction there; also used for the cuts) or fixed values (the correction applied in the experiment, e.g.
  from the GRAPE status). With the GRAPE goal the GRAPE tab's freedom is used, fitted at its points. The
  status after Refresh shows the phase and delay used. Shown as log10(1 − F) or F,
  with worst / mean infidelity and the fraction below 1e-3 / 1e-4 / 1e-5. Not live: Refresh takes the
  current pulse and computes in a separate process (the window does not block and does not interfere with
  a running optimization; a new Refresh replaces a running one). A note appears when the pulse changed.
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
possible (full amplitude) or fixed. Pulse families (2H, 4H, 6H) are generated directly from their phase-ramp
construction. 6H (identity, robust to third order in the offset): three closed 2H loops A B B | A B B, the
second half shifted by 180°; you choose the shape T_A (loop A's helix angle, π < T_A < 4.65, shortest near
4.15), the maximum amplitude and the start phase; loop B and the twist are solved, and the sampled pulse is
polished on the sub-pulse grid (about 1 s), so A1 = A2 = A3 = 0 and the identity hold for the pulse as shown. "calculate pulse sequence live" can be switched off for slow families: then only the curve
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
- Mode "quick": only the gate goal (UR or PP), the duration, the bandwidth (offset range, kHz), the
  maximal amplitude / phase rates and a time budget (minutes). Everything else is the recipe of the
  quick-settings study: rectangular start, broadband grid with 2 points per 1/T, spline controls (one
  per 6 µs) with a small smoothness weight, 22 starts with one screening round (the best third
  continues), the rate limits switched on in a last stage, and in the last 20 % of the time the 3 best
  pulses are re-optimized with the limits. The candidates run in short slices on all cores but one, so
  the run ends on time. Result: the pulse with the smallest worst-case infidelity over the band (401 test
  offsets) that keeps the rate limits (5 % tolerance); the status shows its worst / mean infidelity and
  rates. The polish always gets time for at least 3 slices (the search ends earlier when slices are long,
  e.g. with a B1 range); only pulses better than the start count, and if none is found in time the
  status says so (budget too short). Stop ends early and keeps the best so far. Tuned for UR gates over ±10 kHz at 10 kHz rf (in
  the study: 3–4e-4 worst infidelity for 410–450 µs, 2e-5 for 180° in 600 µs, in 20 min; most of it
  within the first few minutes). Optional B1 range (from–to %, 100 to 100 = off): the pulse should work
  over that amplitude range too. The offsets are then combined with n_B1 = 2 (s_max − s_min) · rf · T + 1
  (at least 3) amplitude scalings (the offset rule applied to the shift (s − 1) · rf of the nutation
  frequency), and the run is scored on 201 offsets × 4 times as many B1 values (the worst value of the area
  decides; the status shows where it is). Each evaluation is ~n_B1 times slower, so a longer budget helps
  (the line below the range shows the point count). The B1 rule is not yet part of the study
  (studies/quick_b1_study.py checks it).
  Freedom: the gate goal's freedom row (none / free phase / linear phase) also applies in quick mode; the
  run is then scored up to the z-rotation fitted at the points, and the status shows the fitted phase and
  delay (not yet part of the study).
  Start pulse: rectangular (default), from a file (.bruker / .csv), 4H or 6H. A file sets the duration (the
  first one of a series); it is resampled to the editor's sub-pulse duration and rf amplitude (same rf field)
  and stretched to other durations. In a series without warm start every duration starts from it. 4H / 6H are
  the robust identity pulses of the pulse families (offset robust to 2nd / 3rd order, 398 / 490 µs at 10 kHz),
  brought to the duration with the same rotations (time stretched, amplitude scaled; shorter than that the
  amplitude is clipped, shown in the info line); GRAPE turns them into the gate.
  With a real start pulse (file, 4H / 6H, warm start) the 22 starts are a mix: the start pulse itself, 10
  copies perturbed only slightly (spread 0.1 instead of 1.0) and 11 random pulses. With the rectangular start
  the recipe is as studied (the perturbed copies are far from it). Start-pulse study
  (studies/start_pulse_study.py): on a wide band (±10 kHz, B1 85–115 %) a 4H start with the mix gained about
  0.2 decades of the worst infidelity over the rectangular quick run; on a narrow band the start did not matter;
  only-near starts ("refine") can get stuck, which the random half avoids.
  Output folder (optional for a single pulse): the result is saved as
  pulse_{µs}_{goal}_{range}_{index}.bruker, e.g. pulse_0550us_UR235-180_off20_B140_2.bruker: duration
  rounded to µs (4 digits); the goal most points have, UR{axis}-{angle} or PP{m0}-{target} with vectors as
  single digits (the smallest integer direction within 0.5°, "m" = minus: 0m10 = −y; angle in whole
  degrees); off{offset range width in kHz}_B1{B1 range width in %} (rounded to 0.1, "p" for the
  decimal point: off7p5); index = the first one that does not
  overwrite an existing file. All numbers rounded half up.
  Series ("series to … , n pulses"): a pulse library, one quick run per duration (evenly spaced, rounded to
  µs); the time is the total, split in proportion to the durations. Warm start (option): from the second
  duration on, the start pulse is the previous result stretched to the new duration (amplitude and
  unwrapped phase interpolated) instead of the chosen start pulse (a real start: the mix above). The pulses go into a new folder named
  like a single pulse with the duration range (pulse_0100-0600us_UR235-180_off20_B140_1/) as
  pulse_0405us.bruker, plus library.json with the settings and the score of every pulse (written after
  each one). Stop ends the series after saving the current pulse (best so far). "Library…" opens the
  library window: worst / mean infidelity vs. duration and a table; a click on a row loads that pulse;
  "Open library…" shows an existing library folder. Mode "expert" shows all the settings below.
- Start pulse: the current pulse of another tab, a rectangular pulse (duration, amplitude, phase) or a
  random one (duration, seed); "Load start pulse". "Result → start" continues from the result.
- Points and goals: the goals hold at a set of points (offset, B1 scaling, weight).
  • offsets [kHz] and B1 [%]: values separated by commas; a:b:n gives n evenly spaced values from a to b
    (e.g. "-10:10:21, 12, 15"). The points are every offset with every B1 value (weight 1).
  • Edit visually…: a separate window with the offset × B1 plane. Click adds a point (optionally also at
    −offset), drag moves the selection, drag on empty space selects in a rectangle or fills it with a grid,
    right click deletes, Shift + wheel changes a point's weight (marker size); exact values in the table.
    Lock to grid snaps to the steps of Settings → Solver (default 0.1 kHz, 0.5 %). Undo / redo, save / load
    (JSON or CSV: offset_kHz, B1_percent, weight). Optional background: log10(1 − F) of the current pulse
    for the gate goal. Apply: a full grid with weights 1 goes back into the fields, anything else is a
    custom selection ("Back to grid" returns to the fields).
  • Goals at every point: the gate goal (order 0; with only that, classical broadband GRAPE) plus the
    ticked offset orders A1/A2/A3 (SCQC) and first-order B1, expanded around every point.
  • Own goals per point (point editor: select points, "Goal for selected…" or double-click the goal in the
    table): a UR rotation or a PP transfer instead of the gate goal of the tab, with presets ("do nothing",
    "excite z → x", "invert z → −z", …). For selective pulses: e.g. excite inside a band, do nothing
    outside. The markers are coloured by goal (white: default). Between the points (area score, background,
    fidelity window) every location is judged by the goal of its nearest point. With own goals for all
    points the gate goal of the tab may be "none".
  • Weights: a point's goal terms count with its weight. All goal terms are divided by √(sum of the
    weights), so the optimizer minimizes a weighted mean over the points: adding points does not shift
    the balance against smoothness, stay close and the rate limits. Lawson starts from the weights.
  • Solved (stop when solved): the gate goal counts as reached when the infidelity is below "'solved' gate
    infidelity below (GRAPE)" (Settings → Solver, default 1e-12 ≈ 1.4e-4° rotation error) at every point;
    the orders when |A_k|/T^k and |B1|/(Ω T) are below "'solved' below".
  • Score over the area (optional, on by default): after the run, the gate infidelity on 201 offsets
    (× 11 B1 values if the points have several) over the range of the points, i.e. also between them;
    the status shows the worst value and where it is, the mean and the fraction below 1e-4.
  • Shortest pulse… (button next to the orders): a small window that searches the shortest pulse doing the
    gate goal exactly and cancelling the offset error through the ticked orders (1st, 1st + 2nd or 1st – 3rd,
    no gaps) and optionally first-order B1, on resonance, at most 100 % rf (the editor's rf amplitude and
    sub-pulse duration; the tab's points and freedom are not used, freedom must be none). Method: GRAPE in
    local mode solves the goals at the start duration (lengthened by 25 % up to 4 times if needed), then the
    duration is shortened in 10 % steps and, after the first failure, by bisection down to one sub-pulse.
    Candidates at every duration: the shortest solution so far squeezed in time (amplitude scaled up, clipped
    at 100 %) and 2 perturbed copies; only if they fail, fresh starts (3 random pulses, CORPSE, 2H, 4H, 6H,
    optionally the tab's start pulse). A candidate that stops improving is given up early. Every accepted
    pulse is checked independently (exact jets: gate error < 1e-10, |A_k|/T^k < 1e-8, B1 < 1e-8) and shown in
    the tab's plots; the shortest one stays as the tab's result ("Result → start", Save… in the window).
    The status gives a simple lower bound (gate only: angle / (360° · rf); with 1st order: 1 / rf, Fenchel's
    theorem). It is a local search: the shortest pulse found, not a proof. Reliability check
    (studies/shortest_reliability.py, 10 kHz, 0.5 µs; Rx(90°) / Rx(180°), orders 1–3, with and without B1,
    3 seeds and start durations each): all 12 cases gave the same minimum 3 times, e.g. Rx(180°) 1st order
    117.0 µs (short CORPSE, 116.7 µs, on the 0.5 µs grid), Rx(90°) 1st order 102.0 µs, 2nd order 201.5 /
    214.5 µs, 3rd order 239.5 / 252.0 µs, with B1 180.0 µs (1st), 214.0 / 228.5 µs (2nd), 289.0 / 288.0 µs
    (3rd); 4–36 min per search on 10 processes. Seed: repeat with another seed and start duration to check a
    result. Time limit: the search ends with the shortest pulse so far (checked between durations); Stop too.
- Gate goal: none, UR (rotation angle about an axis; the whole rotation) or PP (point to point: initial
  magnetization m0 to a target). For PP the robustness terms only count the error that moves m0.
- Freedom (expert mode): the gate goals only have to hold up to a z-rotation after the pulse, fitted by
  the optimizer for every pulse: "free phase" (one angle, absorbed by the phases of later pulses or the
  receiver) or "linear phase" (an angle plus offset × delay, i.e. a free delay after the pulse, correctable
  by a first-order phase correction; |delay| at most "max |delay|", auto = the pulse duration). Often
  allows shorter or better pulses (e.g. broadband excitation). The status shows the fitted phase and delay;
  the area score and the fidelity window use the same fitted rotation.
- The amplitude limit (100 %) is always kept exactly (built into the parametrization).
- Even out points (Lawson reweighting, both modes; off by default): least squares minimizes the sum over
  the points (offsets in local mode, grid points in broadband mode), so the result can be very good at
  some points and worse in between them / at others. If the goals are not reached everywhere, each round
  raises the weights of the points with the largest errors, w ← w·(error / mean)^strength, and
  re-optimizes; this moves towards the smallest worst case (minimax) at the cost of the mean. The round
  with the smallest worst-point error is kept. Strength 0.5 (default) is damped; 1 is classic Lawson
  (faster, may oscillate). "point error spread" in the result = worst / mean point error (1 = even).
- Controls: "every sub-pulse" (classical) or "smooth spline" with M control points per channel (cubic
  B-splines: smooth by construction, few variables, often much faster; needs M ≤ number of sub-pulses).
- Smoothness weight (0 = off; try 1e-4): prefers smooth u_x, u_y (no spikes). It is soft: the optimizer
  first solves with it, then lowers it step by step (×0.1) until the goals are reached exactly, so the
  result stays close to the smooth solution. Too strong a weight can lead into a smooth family without a
  solution: if no candidate reaches the goals with 1e-4 but they do without smoothing, use 1e-5.
- Max. amplitude rate (%/µs) and max. phase rate (kHz = torsion of the error curve / 2π): soft limits that
  cost nothing below the limit (kept to about 1 %). The phase rate counts in proportion to the amplitude
  (a phase jump where the amplitude is ~0 is harmless). Too strict limits make the goals unreachable.
- "Limits only in a last stage" (recommended, on by default): the optimization first runs without the rate
  limits, then a last stage switches them on. Limits from the start can block the way to good solutions,
  which often pass through rough intermediate pulses (example: broadband 180°x ±10 kHz, 450 µs reached a
  worst infidelity of 3e-4 with late limits of 20 %/µs and 50 kHz, but got stuck at 0.67 with the same limits
  from the start). A strong smoothness weight can block the way in the same manner.
- Regularization "stay close to start": adds λ·|u − u_start|² (also lowered step by step).
- Exploration preset: exploit (one start, refined to the end: fastest, finds the nearest solution),
  balanced (8 smooth perturbations, screening 8 → 4 → 2, then 4 kicks) or explore (24 starts, half of
  them random pulses, coarse to fine, screening 24 → 12 → 6 → 3, then 8 kicks with temperature). Changing
  any value gives "custom". The result is the best candidate: goals reached first, then the smoothest (if
  a smoothness weight is set) or the lowest residual. All candidates run in parallel processes (one
  residual line each; dropped candidates turn grey).
- Starts / spread: number of starts and how much starts 2, 3, … differ from the start pulse (0.5 is a
  large perturbation, 0.05 a small one).
- Advanced… (the knobs of exploration vs. exploitation):
  - starts 2, 3, …: perturbations of the start pulse, random pulses, or a mix;
  - noise correlation (µs): smooth noise changes the large-scale shape of the starts and kicks; white noise
    mostly adds spikes that the optimizer then removes;
  - screening rounds (successive halving): every candidate gets "first round" evaluations, the best "kept
    fraction" continue with twice (1 / fraction) the budget, and so on; after the last round the rest run
    to the end. More starts and a smaller first round = more exploration for the same compute;
  - crossover (experimental, genetic): refills a fraction of the dropped places with children of two
    survivors (a blend, or the first part of one and the rest of the other). In tests it did not beat
    screening and kicks alone, so the presets leave it off;
  - kicks: perturb the best result by the kick size, re-optimize, keep it if better. A temperature also
    accepts a worse result sometimes (probability exp(−Δ/T), Δ = increase of log10 of the residual), which
    lets the kicks leave a local minimum. Kicks are skipped once the goals are reached, unless a
    smoothness weight is set (then a smoother solution wins);
  - coarse to fine: first stages on few spline points (e.g. 8, 16), then the chosen controls. Coarse
    stages can only change the large-scale shape, so they explore by construction.
- Start noise: a tiny random change of the first start (seeded). An exactly symmetric start, e.g. a
  constant phase, can otherwise keep the optimizer in a symmetric family of pulses that cannot reach the goal.
The optimizer (Levenberg–Marquardt with exact Jacobians) runs in the background (Optimize / Stop); the plots
show the start (dashed), the best pulse so far and the goal residual vs. evaluations. Stop asks the workers for
their best pulse and waits up to 10 s; if none returns in time (e.g. a busy machine), the best pulse shown so
far is kept. Stopping rules (the
evaluation limit counts per candidate), "solved" thresholds and goal weights are the shared solver settings
("Settings…"; the finite-difference step is not used: the gradients are exact).
If all starts and kicks converge to about the same residual without reaching the goals, the pulse is most
likely too short for them (or the limits too strict): try a longer start pulse.
"""

# the main window's documentation includes the editor description
DOCUMENTATION_STR = DOCUMENTATION_STR.rstrip() + "\n\n\n" + EDITOR_DOCUMENTATION_STR
