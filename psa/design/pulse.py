"""Pulse-sequence data model with CSV and Bruker (JCAMP shape) input/output."""
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np

from psa.getPulseSequence import getPulseSequence


@dataclass
class PulseSequence:
    """Piecewise-constant pulse: N sub-pulses of duration ``dt`` with amplitude (percent of
    ``umax``) and phase (degrees, not wrapped)."""

    amplitude: np.ndarray
    phase: np.ndarray
    dt: float = 0.5e-6
    umax: float = 1e4

    def __post_init__(self):
        self.amplitude = np.asarray(self.amplitude, dtype=np.float64).ravel()
        self.phase = np.asarray(self.phase, dtype=np.float64).ravel()
        if self.amplitude.shape != self.phase.shape:
            raise ValueError("amplitude and phase must have the same length")
        if self.dt <= 0 or self.umax <= 0:
            raise ValueError("dt and umax must be positive")

    @property
    def N(self) -> int:
        return len(self.amplitude)

    @property
    def duration(self) -> float:
        return self.N * self.dt

    def as_array(self) -> np.ndarray:
        """N x 2 array [B1 %, phase deg] as used by the rest of psa."""
        return np.column_stack([self.amplitude, self.phase])

    @classmethod
    def from_array(cls, PS, dt=0.5e-6, umax=1e4) -> "PulseSequence":
        PS = np.asarray(PS, dtype=np.float64).reshape(-1, 2)
        return cls(PS[:, 0].copy(), PS[:, 1].copy(), dt, umax)

    @classmethod
    def load(cls, path, dt=0.5e-6, umax=1e4) -> "PulseSequence":
        """Read a .bruker or .csv file ([B1 %, phase deg] per line). The Bruker header's
        timestep / nominal RF are used when present."""
        PS = getPulseSequence(path)
        if str(path).lower().endswith(".bruker"):
            dt, umax = _bruker_parameters(path, dt, umax)
        return cls.from_array(PS, dt, umax)

    def save_csv(self, path):
        np.savetxt(path, self.as_array(), delimiter=",")

    def save_bruker(self, path, title=None):
        """Bruker JCAMP shape file; phases are wrapped to [0, 360)."""
        title = title or os.path.basename(path)
        phase = np.mod(self.phase, 360.0)
        lines = [
            f"##TITLE= {title}",
            "##JCAMP-DX= 5.00 Bruker JCAMP library",
            "##DATA TYPE= Shape Data",
            "##ORIGIN= psa pulse sequence editor",
            f"##SHAPE_PARAMETERS= Nominal RF (Hz): {self.umax:g}",
            f"##SHAPE_PARAMETERS= timestep: {self.dt:g}",
            f"##MINX= {self.amplitude.min():.6e}" if self.N else "##MINX= 0",
            f"##MAXX= {self.amplitude.max():.6e}" if self.N else "##MAXX= 0",
            f"##MINY= {phase.min():.6e}" if self.N else "##MINY= 0",
            f"##MAXY= {phase.max():.6e}" if self.N else "##MAXY= 0",
            f"##NPOINTS= {self.N}",
            "##XYPOINTS= (XY..XY)",
        ]
        lines += [f"{a:.15f}, {p:.15f}" for a, p in zip(self.amplitude, phase)]
        lines.append("##END=")
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    def save(self, path):
        if str(path).lower().endswith(".bruker"):
            self.save_bruker(path)
        else:
            self.save_csv(path)


def _bruker_parameters(path, dt, umax):
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if not line.startswith("##"):
                continue
            text = line.strip()
            try:
                if "timestep:" in text:
                    dt = float(text.split("timestep:")[1].split()[0])
                elif "Nominal RF (Hz):" in text:
                    umax = float(text.split("Nominal RF (Hz):")[1].split()[0])
            except (ValueError, IndexError):
                pass
    return dt, umax
