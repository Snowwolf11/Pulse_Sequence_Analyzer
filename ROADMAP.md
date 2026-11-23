# Roadmap: ideas for the pulse sequence editor

Collected 2026-11-23. Status: [ ] open, [~] in progress, [x] done.

## A. More kinds of targets
- [x] **A1 Per-point targets (selective pulses).** Every point (offset, B1) may have its own goal (a UR
  rotation or a PP transfer) instead of the gate goal of the GRAPE tab: band-selective excitation /
  inversion, "do nothing" outside a band, solvent suppression, multi-band pulses. Set in the visual point
  editor; default stays the GRAPE tab's gate goal.
- [x] **A2 Goals with freedom.** UR up to a z-rotation; excitation (PP) with a phase that may change
  linearly with the offset (correctable by first-order phase correction, as in BEBOP-type pulses). Often
  allows shorter pulses. Done as a z-rotation after the pulse (free phase, or phase + offset × delay), fitted
  in closed form / by a 1D search at every evaluation (variable projection). Open: free rotation *before*
  the pulse for UR; freedom in the quick mode; "don't care" transition bands for selective pulses.

## B. Hardware and experiment
- [ ] **B3 RF energy / average power** as a soft term or a hard limit (cryoprobes, heating, long sequences).
- [ ] **B4 Constant-amplitude (phase-only) pulses** as an option (adiabatic-like).
- [ ] **B5 Calibration helper:** power level from the hard-90° calibration, integration factor, peak /
  average power of a shape.
- [ ] **B6 Transients and digitization:** amplitude / phase transients, finite amplifier bandwidth, phase /
  amplitude resolution, in the evaluation or in the optimization.

## C. More physics
- [ ] **C7 Relaxation (T1, T2) during the pulse** (long and broadband pulses).
- [ ] **C8 J-coupling / two-spin systems** during the pulse (heteronuclear, long pulses in sequences).
- [ ] **C9 Cooperative pulses (COOP):** pulses designed jointly to work within a phase cycle or a sequence.

## D. Workflow and analysis
- [ ] **D10 Duration scan / own library:** optimize over a range of durations (e.g. the quick recipe), worst /
  mean vs. duration, comparison with a library (e.g. BURBOP).
- [ ] **D11 Phase profile of excitation pulses:** phase of the transverse magnetization vs. offset with a
  linear fit (fidelity window).
- [ ] **D12 Bloch-sphere trajectory animation** over time for chosen offsets / amplitudes.

## Also named for later
- quveros, barq (bigger features, to be specified)

## GUI principle
Keep the main panels simple: details go into small separate setting windows (as the point editor, the
exploration "Advanced…" dialog and the fidelity window) rather than into the main tab.
