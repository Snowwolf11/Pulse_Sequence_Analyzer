"""Quick mode with a B1 range: how many B1 values does the recipe need?

The quick recipe (quick_settings_study.py) was tuned without B1 robustness. For an amplitude range s_lo..s_hi
the quick mode uses n_B1 = round(d (s_hi - s_lo) nu_rf T) + 1 (at least 3) B1 values with d = 2 (the offset
density, by analogy: (s - 1) nu_rf is a frequency like an offset). This study compares d = 1, 2, 3 and, as a
baseline, d = 0 (nominal amplitude only), on three tasks (10 kHz rf, 0.5 µs sub-pulses, rectangular start):

    UR90   UR 90° about x, ±10 kHz, B1 80-120 %, 450 µs   (n_B1 = 3 / 5 / 6 for d = 1 / 2 / 3)
    UR180  UR 180° about x, ±10 kHz, B1 85-115 %, 500 µs  (3 / 4 / 5)
    PPex   PP z -> x, ±10 kHz, B1 80-120 %, 350 µs        (3 / 4 / 5)

Every run is the real QuickRun driver (psa.design.quick) on ONE core with an emulated wall clock of P = 11
cores (as the GUI): per cycle the ready slices run one after another and the clock advances by
max(sum of slice CPU / P, longest slice). Budget: 5 min emulated wall (the GUI default is 10 min; this is a
relative comparison). Scoring on a fixed held-out grid over the area, the same for every variant: 401 offsets x
41 B1 values; main score log10(worst infidelity), also the mean. The incumbent (best pulse within the rate
limits) is also scored at 1.25 / 2.5 / 5 min (anytime curve).

Run:   python studies/quick_b1_study.py --out DIR [--workers 7] [--budget 300] [--smoke]
Report: python studies/quick_b1_study.py --out DIR --report
Interrupted runs resume (results.jsonl is appended; finished jobs are skipped).
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

from psa.design.grape import band_quality, run_task  # noqa: E402
from psa.design.quick import QuickRun, b1_count, quick_config, rectangular_start  # noqa: E402

DT, UMAX, P = 0.5e-6, 1e4, 11
TASKS = {
    "UR90": dict(gate="UR", target_axis=(1, 0, 0), target_angle_deg=90.0, band=(-1e4, 1e4), b1=(0.8, 1.2), T=450e-6),
    "UR180": dict(gate="UR", target_axis=(1, 0, 0), target_angle_deg=180.0, band=(-1e4, 1e4), b1=(0.85, 1.15), T=500e-6),
    "PPex": dict(gate="PP", m0=(0, 0, 1), m_target=(1, 0, 0), band=(-1e4, 1e4), b1=(0.8, 1.2), T=350e-6),
}
DENSITIES = (0.0, 1.0, 2.0, 3.0)
SEEDS = (0, 1)
TRUTH = (401, 41)
CHECKS = (0.25, 0.5, 1.0)                 # fractions of the budget


def jobs(budget):
    out = []
    for task in TASKS:
        for d in DENSITIES:
            for seed in (SEEDS[:1] if d == 0 else SEEDS):
                out.append(dict(key=f"{task}|{d:g}|{seed}|{int(budget)}", task=task, density=d, seed=seed, budget=budget))
    return out


def setup(job):
    t = TASKS[job["task"]]
    kw = {k: t[k] for k in ("gate", "target_axis", "target_angle_deg", "m0", "m_target") if k in t}
    b1 = t["b1"] if job["density"] > 0 else (1.0, 1.0)
    cfg, st = quick_config(t["T"], t["band"], seed=job["seed"], dt=DT, b1_range=b1, umax_hz=UMAX,
                           b1_density=job["density"] or None, **kw)
    return t, cfg, st


def truth(ps, cfg, t):
    q = band_quality(ps, cfg, t["band"][0], t["band"][1], TRUTH[0], t["b1"][0], t["b1"][1], TRUTH[1])
    return {"worst": q["worst infidelity"], "mean": q["mean infidelity"], "worst_at": [q["worst at offset [kHz]"],
            q["worst at B1 [%]"]], "below_1e-3": q["fraction below 0.001"], "below_1e-4": q["fraction below 0.0001"]}


def run_job(job):
    t0 = time.time()
    try:
        t, cfg, st = setup(job)
        wall = [0.0]
        q = QuickRun(rectangular_start(t["T"], DT, UMAX), cfg, st, job["budget"], parallel=P, clock=lambda: wall[0])
        checks = [c * job["budget"] for c in CHECKS]
        snaps, cpu = {}, 0.0
        while not q.done:
            tasks = q.tasks()
            if not tasks:
                break
            times, results = [], []
            for task in tasks:
                c0 = time.process_time()
                results.append(run_task(task))
                times.append(time.process_time() - c0)
            cpu += sum(times)
            before = wall[0]
            wall[0] += max(sum(times) / P, max(times))
            for r in results:
                q.submit(r)
            for c in checks:                                   # incumbent at the checkpoints passed now
                if before < c <= wall[0] and c not in snaps and q.incumbent is not None:
                    snaps[c] = (list(q.incumbent["amp"]), list(q.incumbent["phase"]))
        pulse, s = q.best()
        for c in checks:
            snaps.setdefault(c, (list(pulse.amplitude), list(pulse.phase)))
        from psa.design.pulse import PulseSequence
        anytime = {f"{c:g}": truth(PulseSequence(np.asarray(a), np.asarray(p), DT, UMAX), cfg, t)["worst"]
                   for c, (a, p) in sorted(snaps.items())}
        return dict(job, ok=True, n_b1=int(cfg.bb_scale_n), n_points=len(cfg.ensemble()), truth=truth(pulse, cfg, t),
                    internal_worst=s["band"]["worst infidelity"], within_limits=s["within_limits"], anytime=anytime,
                    slices=s["slices"], wall_emulated=wall[0], cpu_s=cpu, real_s=time.time() - t0,
                    amp=list(map(float, pulse.amplitude)), phase=list(map(float, pulse.phase)))
    except Exception as e:
        return dict(job, ok=False, error=f"{type(e).__name__}: {e}", trace=traceback.format_exc(limit=4))


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
    lines = ["task    d   n_B1  points  seeds  log10 worst (mean of seeds)  mean infid.   worst at 1.25/2.5/5 min"]
    for task in TASKS:
        for d in DENSITIES:
            rs = [r for r in results.values() if r["task"] == task and r["density"] == d]
            if not rs:
                continue
            lw = [np.log10(r["truth"]["worst"]) for r in rs]
            mean = np.mean([r["truth"]["mean"] for r in rs])
            any_ = [np.mean([np.log10(r["anytime"][k]) for r in rs]) for k in sorted(rs[0]["anytime"], key=float)]
            lines.append(f"{task:6s} {d:3g}  {rs[0]['n_b1']:4d}  {rs[0]['n_points']:6d}  {len(rs):5d}  "
                         f"{np.mean(lw):8.2f}  ({', '.join(f'{v:.2f}' for v in lw)})      {mean:9.2e}   "
                         + " / ".join(f"{v:.2f}" for v in any_))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--budget", type=float, default=300.0, help="emulated wall-clock seconds per run")
    ap.add_argument("--smoke", action="store_true", help="tiny budget and one task, to check the pipeline")
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
    budget = 6.0 if a.smoke else a.budget
    todo = [j for j in jobs(budget) if j["key"] not in load(path)]
    if a.smoke:
        todo = [j for j in todo if j["task"] == "PPex" and j["seed"] == 0][:2]
    todo.sort(key=lambda j: -j["density"])                     # the slow ones first
    log = open(os.path.join(a.out, "log.txt"), "a")
    print(f"{time.strftime('%H:%M:%S')} {len(todo)} runs, budget {budget:g} s emulated, {a.workers} workers",
          file=log, flush=True)
    with mp.get_context("spawn").Pool(a.workers) as pool, open(path, "a") as f:
        for r in pool.imap_unordered(run_job, todo):
            f.write(json.dumps(r) + "\n"); f.flush()
            msg = (f"{r['key']}: worst {r['truth']['worst']:.2e}, mean {r['truth']['mean']:.2e}, n_B1 {r['n_b1']}, "
                   f"{r['real_s'] / 60:.1f} min" if r["ok"] else f"{r['key']}: ERROR {r['error']}")
            print(f"{time.strftime('%H:%M:%S')} {msg}", file=log, flush=True)
    text = report(load(path))
    with open(os.path.join(a.out, "report.txt"), "w") as f:
        f.write(text + "\n")
    print(f"{time.strftime('%H:%M:%S')} done\n{text}", file=log, flush=True)


if __name__ == "__main__":
    main()
