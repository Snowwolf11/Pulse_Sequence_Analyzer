#!/usr/bin/env python3
"""
Pulse-sequence parser.

Transforms a domain-specific pulse string into an Nx2 ndarray where:
  - column 0: B1 amplitude in percent
  - column 1: phase in degrees
Each row is one subpulse. Default subpulse duration is 0.5 microseconds.

Supported blocks (end every block with a semicolon ';'):
  - R(B1_percentage, phase, duration);
      duration is optional; defaults to one subpulse (0.5 µs).
  - PR(B1_percentage, phase_start, phase_end, duration);
      duration is optional; phase ramps linearly start->end across subpulses.

Example:
    seq = \"\"\"R(99.1993, 134.201134, 0.5);
    R(100.0,   20, 2.5);
    R(18.0,    -160.0, 10.0);
    PR(99.1, 18.3, 115.4, 20.5);
    R(14.0, 0.0);\"\"\"
    arr = parse_pulse_sequence(seq)
"""

from __future__ import annotations
import re
from typing import List, Tuple
import numpy as np


def parse_pulse_sequence(seq_str: str,
                         subpulse_dt_us: float = 0.5,
                         *,
                         tolerance_us: float = 1e-6) -> np.ndarray:
    """
    Parse a pulse sequence string into an Nx2 ndarray [B1%, phase_deg].

    Parameters
    ----------
    seq_str : str
        Input text containing blocks ending with ';'.
    subpulse_dt_us : float, default 0.5
        Duration of a single subpulse in microseconds.
    tolerance_us : float, default 1e-6
        Allowed deviation when checking duration is a multiple of subpulse_dt_us.

    Returns
    -------
    np.ndarray
        Shape (N, 2). Column 0 = B1 amplitude (%), column 1 = phase (deg).
        If no valid blocks are found, returns an empty array of shape (0, 2).

    Raises
    ------
    ValueError
        On malformed blocks or durations far from multiples of subpulse_dt_us.
    """
    if subpulse_dt_us <= 0:
        raise ValueError("subpulse_dt_us must be positive.")

    # Regex to capture blocks like R(...); or PR(...);
    block_re = re.compile(r'(?is)\b(R|PR)\s*\(\s*([^)]+?)\s*\)\s*;')

    rows: List[Tuple[float, float]] = []

    def _nums_inside(parens: str) -> List[float]:
        # Split by commas; ignore empty entries; allow spaces
        parts = [p.strip() for p in parens.split(',')]
        parts = [p for p in parts if p]
        try:
            return [float(p) for p in parts]
        except ValueError as e:
            raise ValueError(f"Could not parse numeric parameters: '{parens}'") from e

    for m in block_re.finditer(seq_str):
        kind = m.group(1).upper()
        nums = _nums_inside(m.group(2))

        if kind == 'R':
            if len(nums) == 3:
                b1, phase, dur = nums
            elif len(nums) == 2:
                b1, phase = nums
                dur = subpulse_dt_us  # default: one subpulse
            else:
                raise ValueError(f"R block expects 2 or 3 numbers, got {len(nums)}: {nums}")

            if dur < 0:
                raise ValueError(f"Duration must be non-negative in block: {m.group(0)}")

            n = int(round(dur / subpulse_dt_us))
            if abs(n * subpulse_dt_us - dur) > tolerance_us:
                raise ValueError(
                    f"Duration {dur} µs is not a multiple of subpulse_dt_us "
                    f"{subpulse_dt_us} µs within tolerance {tolerance_us} µs."
                )
            if n == 0:
                continue  # zero-length block; skip

            rows.extend((b1, phase) for _ in range(n))

        elif kind == 'PR':
            if len(nums) == 4:
                b1, p_start, p_end, dur = nums
            elif len(nums) == 3:
                b1, p_start, p_end = nums
                dur = subpulse_dt_us  # default: one subpulse
            else:
                raise ValueError(f"PR block expects 3 or 4 numbers, got {len(nums)}: {nums}")

            if dur < 0:
                raise ValueError(f"Duration must be non-negative in block: {m.group(0)}")

            n = int(round(dur / subpulse_dt_us))
            if abs(n * subpulse_dt_us - dur) > tolerance_us:
                raise ValueError(
                    f"Duration {dur} µs is not a multiple of subpulse_dt_us "
                    f"{subpulse_dt_us} µs within tolerance {tolerance_us} µs."
                )
            if n == 0:
                continue

            phases = np.linspace(p_start, p_end, n)  # inclusive endpoints
            rows.extend((b1, float(ph)) for ph in phases)

        else:
            # Should never happen due to regex, but keep for completeness.
            raise ValueError(f"Unknown block type: {kind}")

    if not rows:
        return np.empty((0, 2), dtype=float)

    return np.asarray(rows, dtype=float)
