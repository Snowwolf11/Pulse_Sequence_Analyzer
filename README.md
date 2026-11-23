# Pulse Sequence Analyzer (PSA)

A Python application for analyzing, visualizing, and exporting geometric properties of RF pulse sequences.

The Pulse Sequence Analyzer transforms pulse sequences into spatial trajectories and error curves, allowing investigation of:

- Rotation behavior
- Stability
- Curvature and torsion
- Offset/amplitude robustness
- Quality maps
- Geometric pulse properties

The project includes multiple GUI frontends:

- Tkinter GUI (legacy fallback)
- PySide6 GUI (modern desktop GUI)
- Dash GUI (browser-based interactive GUI)

---

# Quick Start: Step by Step (no experience needed)

This guide assumes you have never done anything like this before. It takes about 15 minutes the first time. You only have to do steps 1–4 once.

## Step 1: Install Python

Python is the free program that runs the Pulse Sequence Analyzer. The analyzer needs **Python 3.13**.
(Do not use the big yellow "Download Python" button on python.org: it gives a newer version that does not work with the analyzer yet.)

1. Click the link for your computer. The download starts right away:
   - **Mac** (macOS 12 or newer): [python-3.13.17-macos11.pkg](https://www.python.org/ftp/python/3.13.17/python-3.13.17-macos11.pkg)
   - **Windows** (64-bit): [python-3.13.17-amd64.exe](https://www.python.org/ftp/python/3.13.17/python-3.13.17-amd64.exe)
2. Open the downloaded file (it is in your **Downloads** folder) and click through the installer.
   - **Windows only, important:** on the very first screen, tick the box **"Add python.exe to PATH"** before clicking "Install Now".

## Step 2: Download the Pulse Sequence Analyzer

1. On this GitHub page, click the green **"Code"** button near the top.
2. Click **"Download ZIP"**.
3. Open your **Downloads** folder and double-click the ZIP file. A folder called **`Pulse_Sequence_Analyzer-main`** appears.
4. Optional: move this folder somewhere you will find it again, e.g. your Desktop or Documents.

## Step 3: Open the "Terminal" (a window where you type commands)

- **Mac:** press `Cmd` + `Space`, type **Terminal**, press `Enter`.
- **Windows:** click the Start button, type **cmd**, press `Enter` ("Command Prompt").

A window with text appears. Do not worry: you only copy and paste a few lines into it.
After each line, press `Enter` and wait until the window shows a new empty line.

## Step 4: Set up the analyzer (only the first time)

1. In the Terminal, type `cd` followed by **one space** (do not press `Enter` yet).
   On **Windows**, type `cd /d` followed by one space instead.
2. Drag the folder `Pulse_Sequence_Analyzer-main` from Finder / Explorer **into the Terminal window** and let go.
   Its location now appears after `cd `. Press `Enter`.
3. Copy and paste the following lines, one at a time, pressing `Enter` after each:

   **Mac:**
   ```bash
   python3.13 -m venv .venv
   source .venv/bin/activate
   pip install -e .
   ```

   **Windows:**
   ```bat
   py -3.13 -m venv .venv
   .venv\Scripts\activate
   pip install -e .
   ```

   The last line downloads everything the analyzer needs. This takes a few minutes and prints a lot of text. That is normal.
   It is finished when you see a line starting with `Successfully installed` and a new empty line.
   A message like `[notice] A new release of pip is available` may follow. You can ignore it.

## Step 5: Start the analyzer

Type this and press `Enter`:

```bash
python main.py --gui pyside6
```

After a few seconds the Pulse Sequence Analyzer window opens.
Keep the Terminal window open while you use the analyzer; closing it also closes the analyzer.

## Every next time

You do not need to repeat steps 1–4. Just:

1. Open the Terminal (Step 3).
2. Type `cd ` (Windows: `cd /d `), drag the `Pulse_Sequence_Analyzer-main` folder into the window, press `Enter`.
3. Paste these two lines (pressing `Enter` after each):
   - **Mac:** `source .venv/bin/activate` and then `python main.py --gui pyside6`
   - **Windows:** `.venv\Scripts\activate` and then `python main.py --gui pyside6`

## If something goes wrong

| What you see | What to do |
|---|---|
| `command not found: python3.13` (Mac), `'py' is not recognized` or `No suitable Python runtime found` (Windows) | Python is not installed correctly. Repeat Step 1. On Windows, make sure you ticked "Add python.exe to PATH". |
| A Mac pop-up asks to install "command line developer tools" | Click "Cancel" and install Python from python.org as in Step 1 instead. |
| `requires a different Python` or `No matching distribution found for PySide6` | The setup used a Python version newer than 3.13. Install Python 3.13 (Step 1), delete the folder `.venv` inside `Pulse_Sequence_Analyzer-main` (it may be hidden: on Mac press `Cmd`+`Shift`+`.` in Finder to show it), and repeat Step 4. |
| `No such file or directory: main.py` | The Terminal is not "in" the right folder. Repeat Step 4.1–4.2 (the `cd` + drag step). |
| The analyzer window does not open, no error | Try the simpler window instead: `python main.py` |

---

# Features

## Pulse Sequence Analysis

Analyze pulse sequences and compute:

- Error curves
- Magnetization trajectories
- Curvature
- Torsion
- Integrated geometric quantities
- Rotation axis and angle
- Pulse statistics

---

## Interactive Visualization

Interactive plotting of:

- 3D error curves
- 3D trajectories
- Curvature vs arc length
- Torsion vs arc length
- Stability surfaces
- Quality images

---

## Stability Analysis

Compute robustness against:

- Amplitude scaling
- Frequency offset
- Rotation angle deviations

Supports:

- PP mode
- UR mode
- Angle mode

---

## Quality Image Generation

Generate and export quality images for entire directories of pulse sequences.

Supports:

- Offset scans
- Amplitude scans
- SS calculations
- UR calculations

---

## Pulse Sequence Editor (PySide6 GUI)

Create pulse sequences visually, from curves, or from parametrized curve families, with live evaluation of
the gate and its robustness, and a solver that adjusts the design for robustness and an optional target gate.
See [Pulse Sequence Editor](#pulse-sequence-editor) below.

---

## Fast Calculations

All calculations use a vectorized NumPy backend (no compiled extensions needed).

---

# Project Structure

```text
psa/
├── main.py                  # entry point: python main.py [--gui tkinter|pyside6|dash]
├── psa/
│   ├── GUI.py               # Tkinter GUI
│   ├── GUI_pyside6.py       # PySide6 GUI (recommended)
│   ├── GUI_dash.py          # Dash browser GUI
│   ├── app_backend.py       # GUI-independent analysis functions
│   ├── rotations.py         # vectorized NumPy rotation backend
│   ├── createCurve.py, ...  # curve, stability and quality-image calculations
│   ├── design/              # editor back end: pulse <-> curve, metrics, families, builder, solver, GRAPE
│   ├── editor_gui/          # editor front end (PySide6 + pyqtgraph)
│   └── tests/               # regression, physics, editor back-end and GUI tests
├── pyproject.toml
├── requirements.txt
└── environment.yml
```

---

# Installation

Only Python 3.8 – 3.13 is needed; there is nothing to compile.

```bash
git clone https://github.com/Snowwolf11/Pulse_Sequence_Analyzer.git
cd Pulse_Sequence_Analyzer
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .
```

Or with conda:

```bash
conda env create -f environment.yml
conda activate psa
```

---

# Starting the Application

```bash
python main.py --gui pyside6       # recommended desktop GUI
python main.py                     # Tkinter GUI (fallback)
python main.py --gui dash          # browser GUI, then open http://127.0.0.1:8050
```

---

# Running the Tests

```bash
pip install -e ".[test]"
python -m pytest psa/tests
```

`test_physics.py` checks the rotation backend against exact SU(2) propagation;
`test_regression.py` checks all calculations against reference results of the original implementation.

---

# Using the Application

## 1. Load a Pulse Sequence

You can:

- Enter a pulse sequence file
- Browse for a file
- Create/edit sequences using the editor

Supported formats:

- .csv
- .bruker
- generated pulse sequences

---

## 2. Adjust Parameters

Main parameters:

| Parameter | Description |
|---|---|
| Pulse Amplitude | Maximum RF amplitude |
| Vector Length | Magnetization vector length |
| Time/Pulse | Pulse duration |
| Offset | Frequency offset |
| Scaling | RF scaling |
| InpoFact | Interpolation factor |
| x_Expand | Visual x-axis expansion |

---

## 3. Refresh / Recalculate

Press:

- Refresh button
- ENTER key

to recompute all plots and analysis.

---

## 4. Interactive Sliders

Use sliders for:

- Amplitude scaling
- Offset variation

Playback buttons automatically sweep parameters to visualize curve evolution dynamically.

---

## 5. Stability Calculations

Open the Stability tab and configure:

- Amplitude range
- Offset range
- Stability mode

Then click:

```text
Calculate Stability
```

---

## 6. Quality Images

Configure:

- Input directory
- Output directory
- Resolution
- Scan range
- Calculation type

Then calculate and export quality images.

---

# Pulse Sequence Editor

Start the PySide6 interface (`python main.py --gui pyside6`) and click **Edit**. The editor has tabs on the
left and a live panel on the right; its **Help** button explains everything in detail.

| Tab | What you do |
|---|---|
| Text | Write the pulse with the `R(...)` / `PR(...)` syntax (see below). |
| Waveform | Draw amplitude and phase over time with draggable keyframes (linear or hold between them). |
| Curve families | Pick a family (helix, Lissajous knot, torus knot, Clelia curve, custom formula, 2H, 4H, 6H) and change its parameters live. The 6H family is an identity pulse robust to third order in the offset. |
| Curve builder | Compose the error curve from lines, arcs, helices and Bézier segments with a twist (phase jump) at each junction; edit in 3D and in three projections. |
| Import curve | Load curve points (x, y, z) or a curvature/torsion table from a CSV / text file. |
| Optimize (GRAPE) | Optimize every sub-pulse of a start pulse with exact gradients: local robustness (detuning orders and B1 at chosen offsets, as in SCQC) or classical broadband GRAPE, with a UR or PP gate goal; a quick mode (gate, duration, bandwidth, time budget) uses the recipe of the settings study. |

The panel shows the gate (rotation angle and axis), the robustness against detuning to first, second and
third order (end point, area and third-order term of the error curve), the first-order robustness against
amplitude (B1) scaling, the error curve in 3D and the fidelity vs. offset. Pulses are saved as CSV or Bruker
shape files and can be sent directly to the analyzer.

**Solver.** In the Waveform, Curve families and Curve builder tabs a solver adjusts free variables
(keyframe values, family parameters, or builder segment parameters marked "free") so that the chosen goals
are met: cancel the detuning error to 1st, 2nd and/or 3rd order, B1 robustness, and optionally a target
gate (angle about an axis). It runs in the background and can be stopped at any time (the best result so
far is kept). Evaluation limit, tolerances and goal weights are set in Settings → Solver.

**GRAPE.** The Optimize (GRAPE) tab changes every sub-pulse (Cartesian controls) instead of a few
parameters. The goals hold at a set of points (offset, B1 scaling, weight): an offset × B1 grid typed as
lists and ranges (`-10:10:21`), or a custom selection drawn in a visual point editor (weights, lock to grid,
save / load, the infidelity of the current pulse as background). At every point the gate goal is required
(order 0, i.e. classical broadband GRAPE) and optionally the Taylor orders of the detuning error and the
first-order B1 error (SCQC). The
gate goal is none, UR (rotation) or PP (point to point); the 100 % amplitude limit is always kept.
Optional Lawson reweighting evens out the errors across the offsets / grid points (smaller worst case).
Against spiky results: a smoothness weight, soft limits of the amplitude and phase rates, or smooth spline
controls (optional; every sub-pulse stays the default). Exploration vs. exploitation is set with a
preset (exploit / balanced / explore) or in detail: starts (perturbed, smooth or white, or random),
screening rounds that keep the best candidates and give them more budget, kicks of the best result
(optionally with a temperature), coarse-to-fine stages and an experimental crossover. All candidates run in
parallel and the best result is kept. "Stay close to start" regularization is optional. It uses the same solver
settings.
A *quick* mode needs only the gate, duration, bandwidth, rate limits and a time budget; it runs the recipe
found in a settings study (`studies/quick_settings_study.py`: broadband grid, spline controls, 22 starts with
screening, rate limits in a last stage and a final polish of the 3 best) on all cores and returns the pulse
with the smallest worst-case infidelity over the band that keeps the limits.

Conventions: on resonance, unit-speed error curve whose length is the pulse duration, curvature = Rabi
frequency, torsion = phase change rate. Curves are converted into pulses exactly (the pulse reproduces the
curve's tangent at every sub-pulse boundary).

---

# Pulse Sequence Syntax

Example:

```text
R(100, 0, 50); R(100, 90, 50);  PR(100, 0, 360, 90);
```

## Commands

### Rectangular Pulse

```text
R(amplitude [%], phase [°], duration [µs]);
```

### Phase Ramp Pulse

```text
PR(amplitude [%], start phase [°], end phase [°], duration [µs]);
```

---

# Recommended GUI

## Recommended: PySide6 GUI

Advantages:

- Modern dark theme
- Better responsiveness
- Better scalability
- Better interactivity
- Easier future extension

---

# Dependencies

Main dependencies:

- numpy
- scipy
- matplotlib
- plotly
- dash
- PySide6
- pyqtgraph
- mplcursors
- mpld3

---

# Notes

- Tkinter GUI remains available as fallback.
- PySide6 is the primary desktop GUI moving forward.
- Dash GUI enables browser-based visualization and remote usage.

---

# Author

Leon Eppinger