"""Shortest robust pulse search from the command line (psa.design.shortest).

Example:  python studies/shortest_pulse.py --angle 180 --order 1 --start 200 --workers 4 --out x180_o1
          python studies/shortest_pulse.py --angle 90 --axis 0 1 0 --order 2 --b1 --start 400
Writes OUT.bruker (the pulse) and OUT.json (verification + descent history) if --out is given.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import multiprocessing as mp
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from psa.design.shortest import ShortestSpec, search  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", default="UR", choices=["UR", "PP"])
    ap.add_argument("--axis", type=float, nargs=3, default=[1, 0, 0])
    ap.add_argument("--angle", type=float, default=180.0)
    ap.add_argument("--m0", type=float, nargs=3, default=[0, 0, 1])
    ap.add_argument("--mt", type=float, nargs=3, default=[1, 0, 0])
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--b1", action="store_true")
    ap.add_argument("--umax", type=float, default=1e4)
    ap.add_argument("--dt", type=float, default=0.5, help="µs")
    ap.add_argument("--start", type=float, default=300.0, help="start duration (µs)")
    ap.add_argument("--random", type=int, default=3, help="random starts per duration")
    ap.add_argument("--nfev", type=int, default=3000)
    ap.add_argument("--time", type=float, default=900.0, help="time limit (s)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out")
    a = ap.parse_args()
    spec = ShortestSpec(gate=a.gate, axis=a.axis, angle_deg=a.angle, m0=a.m0, mt=a.mt, order=a.order, b1=a.b1,
                        dt=a.dt * 1e-6, umax=a.umax)
    kw = dict(n_random=a.random, max_nfev=a.nfev, time_limit_s=a.time, log=lambda s: print(s, flush=True))
    if a.workers > 1:
        with mp.get_context("spawn").Pool(a.workers) as pool:
            res = search(spec, a.start, mapper=pool.map, **kw)
    else:
        res = search(spec, a.start, **kw)
    v = res["verification"]
    print(f"\n{res['message']}; {res['elapsed_s']:.0f} s")
    if v is None:
        return
    print(f"shortest found: {v['duration [µs]']:.1f} µs (N = {v['N']}), lower bound {res['lower bound [µs]']:.1f} µs")
    print("verification: " + ", ".join(f"{k} {val:.2e}" if isinstance(val, float) else f"{k} {val}" for k, val in v.items()))
    if a.out:
        res["pulse"].save_bruker(a.out + ".bruker", title=f"shortest {a.gate} order {a.order}{' B1' if a.b1 else ''}")
        with open(a.out + ".json", "w") as f:
            json.dump({"spec": spec.__dict__, "verification": v, "history": res["history"],
                       "lower bound [µs]": res["lower bound [µs]"], "message": res["message"]}, f, indent=1)


if __name__ == "__main__":
    main()
