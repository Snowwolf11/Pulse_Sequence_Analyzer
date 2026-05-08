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
Define a pulse sequence as a list of rectangular pulses and phase ramps, using the format shown in the editor when it opens. Choose a target directory; the defined sequence will be saved there. Use **Plot** to visualize the sequence in the Curve Window. Calculations use the pulse parameters from the Main Window.

Menu / Settings
- Calculation method: **Rotation Matrix** or **Helix**.
- Initial condition: select the **error direction**.
- Backend / language:
  • **Rust** — much faster but requires additional setup.
  • **Python** — runs out of the box but is slower.

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
- Note: with the Python backend, this computation can take several minutes.

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
Note: with default settings and ~150 sequences, the Python backend may require several minutes or more.

Tips
- Press **Enter** in any numeric field to recompute quickly.
"""