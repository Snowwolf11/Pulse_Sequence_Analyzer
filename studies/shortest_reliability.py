"""Is the shortest robust pulse search (psa.design.shortest) reliable? Every case is searched 3 times with
different seeds and start durations; a case counts as reproducible when the minima agree within 1 µs.

Cases: Rx(90°), Rx(180°) x order 1, 2, 3 x without / with first-order B1 (12 cases, 36 searches), 10 kHz,
0.5 µs. Known answer: Rx(180°) order 1 without B1 is short CORPSE, 116.7 µs -> 117.0 µs on the grid.
One search per process (sequential candidates inside), time limit per search.

Run:    python studies/shortest_reliability.py --out DIR [--workers 10] [--time 3600]
Report: python studies/shortest_reliability.py --out DIR --report
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import multiprocessing as mp
import sys
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from psa.design.shortest import ShortestSpec, search  # noqa: E402

ANGLES = (90.0, 180.0)
ORDERS = (1, 2, 3)
B1 = (False, True)
REPS = ((0, 150.0), (1, 250.0), (2, 200.0))          # (seed, start duration µs)
AGREE_US = 1.0


def jobs():
    return [dict(key=f"x{a:g}|o{o}|{'b1' if b else 'nob1'}|s{s}", angle=a, order=o, b1=b, seed=s, start=T)
            for a in ANGLES for o in ORDERS for b in B1 for s, T in REPS]


def run_job(job, time_limit):
    t0 = time.time()
    try:
        spec = ShortestSpec(angle_deg=job["angle"], order=job["order"], b1=job["b1"])
        res = search(spec, job["start"], seed=job["seed"], time_limit_s=time_limit)
        v = res["verification"]
        out = dict(job, ok=True, message=res["message"], history=res["history"], real_s=time.time() - t0,
                   verification=v, duration_us=None if v is None else v["duration [µs]"])
        if res["pulse"] is not None:
            out["amp"], out["phase"] = list(map(float, res["pulse"].amplitude)), list(map(float, res["pulse"].phase))
        return out
    except Exception as e:
        return dict(job, ok=False, error=f"{type(e).__name__}: {e}", trace=traceback.format_exc(limit=4))


def _run(args):
    return run_job(*args)


def load(path):
    out = {}
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    if r.get("ok"):
                        out[r["key"]] = r
                except json.JSONDecodeError:
                    pass
    return out


def report(results):
    lines = ["case               minima of the repeats [µs]      spread  reproducible  time [min]  notes"]
    for a in ANGLES:
        for o in ORDERS:
            for b in B1:
                rs = sorted((r for r in results.values() if r["angle"] == a and r["order"] == o and r["b1"] == b),
                            key=lambda r: r["seed"])
                if not rs:
                    continue
                d = [r["duration_us"] for r in rs]
                found = [x for x in d if x is not None]
                spread = (max(found) - min(found)) if found else float("nan")
                notes = sorted({r["message"] for r in rs if "time limit" in r["message"] or r["duration_us"] is None})
                lines.append(f"Rx({a:g}) o{o} {'B1  ' if b else '    '}  "
                             f"{', '.join('-' if x is None else f'{x:6.1f}' for x in d):30s}  {spread:6.1f}  "
                             f"{'yes' if found and len(found) == len(d) and spread <= AGREE_US else 'NO ':12s}  "
                             f"{np.mean([r['real_s'] for r in rs]) / 60:9.1f}  {'; '.join(notes)}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--time", type=float, default=3600.0, help="time limit per search (s)")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "results.jsonl")
    if a.report:
        text = report(load(path))
        print(text)
        with open(os.path.join(a.out, "report.txt"), "w") as f:
            f.write(text + "\n")
        return
    todo = [j for j in jobs() if j["key"] not in load(path)]
    todo.sort(key=lambda j: (-j["order"], -j["b1"]))          # the slow ones first
    log = open(os.path.join(a.out, "log.txt"), "a")
    print(f"{time.strftime('%H:%M:%S')} {len(todo)} searches, {a.workers} workers, {a.time:g} s each at most",
          file=log, flush=True)
    with mp.get_context("spawn").Pool(a.workers) as pool, open(path, "a") as f:
        for r in pool.imap_unordered(_run, [(j, a.time) for j in todo]):
            f.write(json.dumps(r) + "\n"); f.flush()
            msg = (f"{r['key']}: {r['duration_us']} µs ({r['message']}), {r['real_s'] / 60:.1f} min" if r["ok"]
                   else f"{r['key']}: ERROR {r['error']}")
            print(f"{time.strftime('%H:%M:%S')} {msg}", file=log, flush=True)
    text = report(load(path))
    with open(os.path.join(a.out, "report.txt"), "w") as f:
        f.write(text + "\n")
    print(f"{time.strftime('%H:%M:%S')} done\n{text}", file=log, flush=True)


if __name__ == "__main__":
    main()
