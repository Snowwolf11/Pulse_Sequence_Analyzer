"""Keyframe model of the visual amplitude / phase editor.

Each channel (amplitude in %, phase in degrees, phase not wrapped) is a list of keyframes at
sub-pulse indices 0 ... N-1. Between two keyframes the values are interpolated per sub-pulse,
either linearly or held constant ("hold"); the mode belongs to the segment that starts at a
keyframe. Before the first / after the last keyframe the value is held.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np

from psa.design.pulse import PulseSequence

LINEAR, HOLD = "linear", "hold"


@dataclass
class Keyframe:
    index: int
    value: float
    mode: str = LINEAR          # interpolation of the segment starting here


@dataclass
class Channel:
    keyframes: List[Keyframe] = field(default_factory=list)
    minimum: float = -np.inf
    maximum: float = np.inf

    def sort(self):
        self.keyframes.sort(key=lambda k: k.index)

    def clamp(self, value):
        return float(min(self.maximum, max(self.minimum, value)))

    def find(self, index):
        for i, k in enumerate(self.keyframes):
            if k.index == index:
                return i
        return None

    def set(self, index, value, mode=None):
        """Add or replace the keyframe at ``index``."""
        value = self.clamp(value)
        i = self.find(index)
        if i is None:
            self.keyframes.append(Keyframe(int(index), value, mode or LINEAR))
            self.sort()
        else:
            self.keyframes[i].value = value
            if mode:
                self.keyframes[i].mode = mode

    def move(self, position, index, value, N):
        """Move keyframe number ``position`` to (index, value); the index is snapped to 0..N-1 and may
        not pass a neighbouring keyframe. Returns the new position in the list."""
        k = self.keyframes[position]
        lo = self.keyframes[position - 1].index + 1 if position > 0 else 0
        hi = self.keyframes[position + 1].index - 1 if position + 1 < len(self.keyframes) else N - 1
        k.index = int(min(max(int(round(index)), lo), max(lo, hi)))
        k.value = self.clamp(value)
        return position

    def remove(self, position):
        if len(self.keyframes) > 1:
            del self.keyframes[position]

    def values(self, N) -> np.ndarray:
        out = np.empty(N)
        kf = [k for k in self.keyframes if 0 <= k.index < N]
        if not kf:
            out[:] = self.clamp(0.0)
            return out
        out[:kf[0].index] = kf[0].value
        for a, b in zip(kf[:-1], kf[1:]):
            i = np.arange(a.index, b.index)
            out[a.index:b.index] = a.value if a.mode == HOLD else a.value + (b.value - a.value) * (i - a.index) / (b.index - a.index)
        out[kf[-1].index:] = kf[-1].value
        return out

    def truncate(self, N):
        self.keyframes = [k for k in self.keyframes if k.index < N] or [Keyframe(0, self.clamp(0.0))]

    def to_dict(self):
        return [[k.index, k.value, k.mode] for k in self.keyframes]

    @classmethod
    def from_list(cls, data, minimum=-np.inf, maximum=np.inf):
        return cls([Keyframe(int(i), float(v), str(m)) for i, v, m in data], minimum, maximum)


def simplify(values, tolerance):
    """Keyframe indices (Ramer-Douglas-Peucker on index/value) such that linear interpolation between
    them reproduces ``values`` within ``tolerance``. tolerance = 0 keeps every non-linear point."""
    v = np.asarray(values, dtype=np.float64)
    n = len(v)
    if n <= 2:
        return list(range(n))
    keep = np.zeros(n, dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        i = np.arange(a + 1, b)
        dev = np.abs(v[i] - (v[a] + (v[b] - v[a]) * (i - a) / (b - a)))
        j = int(np.argmax(dev))
        if dev[j] > tolerance + 1e-12:
            keep[a + 1 + j] = True
            stack += [(a, a + 1 + j), (a + 1 + j, b)]
    return list(np.flatnonzero(keep))


class Waveform:
    """Amplitude and phase keyframe channels plus metadata; produces a PulseSequence."""

    def __init__(self, N=200, dt=0.5e-6, umax=1e4, amplitude=100.0, phase=0.0):
        self.N, self.dt, self.umax = int(N), float(dt), float(umax)
        self.amplitude = Channel([Keyframe(0, amplitude)], 0.0, 100.0)
        self.phase = Channel([Keyframe(0, phase)])

    @property
    def duration(self):
        return self.N * self.dt

    def set_duration(self, duration):
        """Change the number of sub-pulses (keyframes beyond the new end are removed)."""
        self.N = max(1, int(round(duration / self.dt)))
        self.amplitude.truncate(self.N)
        self.phase.truncate(self.N)

    def set_dt(self, dt):
        """Change the sub-pulse duration keeping the total duration and the keyframe times
        (indices are rescaled; keyframes that collide on one sub-pulse are merged)."""
        factor = self.dt / float(dt)
        self.dt = float(dt)
        self.N = max(1, int(round(self.N * factor)))
        for ch in (self.amplitude, self.phase):
            merged = {}
            for k in ch.keyframes:
                k.index = min(self.N - 1, int(round(k.index * factor)))
                merged[k.index] = k
            ch.keyframes = [merged[i] for i in sorted(merged)]

    def to_pulse(self) -> PulseSequence:
        return PulseSequence(self.amplitude.values(self.N), self.phase.values(self.N), self.dt, self.umax)

    @classmethod
    def from_pulse(cls, ps: PulseSequence, amplitude_tolerance=0.0, phase_tolerance=0.0) -> "Waveform":
        """Keyframes that reproduce ``ps`` (exactly for tolerance 0). The phase is unwrapped first."""
        w = cls(ps.N, ps.dt, ps.umax)
        amp = np.clip(ps.amplitude, 0.0, 100.0)
        phase = np.degrees(np.unwrap(np.radians(ps.phase)))
        w.amplitude = Channel([Keyframe(int(i), float(amp[i])) for i in simplify(amp, amplitude_tolerance)], 0.0, 100.0)
        w.phase = Channel([Keyframe(int(i), float(phase[i])) for i in simplify(phase, phase_tolerance)])
        return w

    def to_dict(self):
        return {"type": "waveform", "N": self.N, "dt": self.dt, "umax": self.umax,
                "amplitude": self.amplitude.to_dict(), "phase": self.phase.to_dict()}

    @classmethod
    def from_dict(cls, d):
        w = cls(d["N"], d["dt"], d["umax"])
        w.amplitude = Channel.from_list(d["amplitude"], 0.0, 100.0)
        w.phase = Channel.from_list(d["phase"])
        return w
