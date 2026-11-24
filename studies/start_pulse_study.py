"""Does the start pulse matter for the quick mode? Rectangular vs BARQ (qurveros, gate-matched) vs 4H / 6H (identity).

Start pulses (10 kHz rf, 0.5 µs sub-pulses), brought to the task duration T gate-preservingly: time stretched,
amplitude scaled by T_natural / T (clipped at 100 %):
    rect   rectangular, 100 %, duration T (the quick-mode default)
    barq   BARQ pulse for the task gate (B1 robust to 1st order, offset by closure, on resonance); one BARQ
           optimization per seed, files UR{angle}_s{seed}.npz in --starts (made with qurveros, see barq_make.py)
    4H     4H 0° UR (τ = 99.5 µs, 398 µs), identity gate: robust shape, wrong gate
    6H     6H 0° UR (T_A = 4.15, 490 µs), identity gate

Modes:
    quick    the quick recipe as is: 22 starts = the start pulse + 10 strongly perturbed copies (spread 1.0)
             + 11 random pulses, screening, best third continues; the start pulse is 1 of 22 candidates
    mix      the quick recipe, but the 10 perturbed copies near the start pulse (perturbation 0.1): 11 near + 11 random
    refine   6 starts, all near the start pulse (perturbation 0.1, "perturb" pool), otherwise the same:
             is the start pulse a good GRAPE start?

Tasks (quick recipe with B1 density 2, as the GUI):
    UR180w   UR 180° x, ±10 kHz, B1 85-115 %, 500 µs   (wide: GRAPE's home ground)
    UR180n   UR 180° x, ±3 kHz,  B1 90-110 %, 500 µs   (narrow: where BARQ / helix pulses are robust already)

Driver, emulated clock, truth grid (401 x 41) and checkpoints as quick_b1_study.py.
Run:    python studies/start_pulse_study.py --out DIR --starts DIR [--workers 7] [--budget 300]
Report: python studies/start_pulse_study.py --out DIR --report
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

from psa.design.families import family_by_name  # noqa: E402
from psa.design.grape import band_quality, run_task  # noqa: E402
from psa.design.pulse import PulseSequence  # noqa: E402
from psa.design.quick import QuickRun, quick_config, rectangular_start  # noqa: E402

DT, UMAX, P = 0.5e-6, 1e4, 11
TASKS = {
    "UR180w": dict(angle=180.0, band=(-1e4, 1e4), b1=(0.85, 1.15), T=500e-6),
    "UR180n": dict(angle=180.0, band=(-3e3, 3e3), b1=(0.9, 1.1), T=500e-6),
}
STARTS = ("rect", "barq", "4H", "6H")
MODES = ("quick", "mix", "refine")
MODE_SETTINGS = {"mix": dict(start_pool="mix", multistart_noise=0.1),          # 22 starts: 11 near the start, 11 random
                 "refine": dict(n_starts=6, start_pool="perturb", multistart_noise=0.1)}
SEEDS = (0, 1)
TRUTH = (401, 41)
CHECKS = (0.25, 0.5, 1.0)


def fit_duration(ps, T):
    """Gate-preserving change of duration: time stretched, amplitude x T_old / T (clipped at 100 %)."""
    n = max(4, int(round(T / DT)))
    x_old, x_new = (np.arange(ps.N) + 0.5) / ps.N, (np.arange(n) + 0.5) / n
    amp = np.interp(x_new, x_old, np.asarray(ps.amplitude, float)) * ps.duration / (n * DT)
    phase = np.degrees(np.interp(x_new, x_old, np.unwrap(np.radians(np.asarray(ps.phase, float)))))
    return PulseSequence(np.clip(amp, 0, 100), (phase + 180) % 360 - 180, DT, UMAX)


def start_pulse(job, starts_dir):
    t = TASKS[job["task"]]
    if job["start"] == "rect":
        return rectangular_start(t["T"], DT, UMAX)
    if job["start"] == "barq":
        d = np.load(os.path.join(starts_dir, f"UR{t['angle']:g}_s{job['seed']}.npz"))
        ps = PulseSequence(d["amp"], d["phase"], DT, UMAX)
    elif job["start"] == "4H":
        ps = family_by_name("4H (0° UR, 2nd-order robust)").pulse({"tau": 99.5}, DT, UMAX)
    else:
        ps = family_by_name("6H (0° UR, 3rd-order robust)").pulse({"TA": 4.15}, DT, UMAX)
    return fit_duration(ps, t["T"])


def jobs(budget):
    return [dict(key=f"{task}|{m}|{s}|{seed}|{int(budget)}", task=task, mode=m, start=s, seed=seed, budget=budget)
            for task in TASKS for m in MODES for s in STARTS for seed in SEEDS]


def truth(ps, cfg, t):
    q = band_quality(ps, cfg, t["band"][0], t["band"][1], TRUTH[0], t["b1"][0], t["b1"][1], TRUTH[1])
    return {"worst": q["worst infidelity"], "mean": q["mean infidelity"]}


def run_job(job, starts_dir):
    t0 = time.time()
    try:
        t = TASKS[job["task"]]
        cfg, st = quick_config(t["T"], t["band"], seed=job["seed"], dt=DT, b1_range=t["b1"], umax_hz=UMAX,
                               gate="UR", target_axis=(1, 0, 0), target_angle_deg=t["angle"])
        if job["mode"] in MODE_SETTINGS:
            for k, v in MODE_SETTINGS[job["mode"]].items():
                setattr(cfg, k, v)
            cfg.check()
        start = start_pulse(job, starts_dir)
        start_truth = truth(start, cfg, t)
        wall = [0.0]
        q = QuickRun(start, cfg, st, job["budget"], parallel=P, clock=lambda: wall[0])
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
            for c in checks:
                if before < c <= wall[0] and c not in snaps and q.incumbent is not None:
                    snaps[c] = (list(q.incumbent["amp"]), list(q.incumbent["phase"]))
        pulse, s = q.best()
        for c in checks:
            snaps.setdefault(c, (list(pulse.amplitude), list(pulse.phase)))
        anytime = {f"{c:g}": truth(PulseSequence(np.asarray(a), np.asarray(p), DT, UMAX), cfg, t)["worst"]
                   for c, (a, p) in sorted(snaps.items())}
        return dict(job, ok=True, start_truth=start_truth, truth=truth(pulse, cfg, t), anytime=anytime,
                    optimized=bool(s.get("optimized", True)), within_limits=s["within_limits"],
                    wall_emulated=wall[0], cpu_s=cpu, real_s=time.time() - t0,
                    amp=list(map(float, pulse.amplitude)), phase=list(map(float, pulse.phase)))
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
    lines = ["task     mode    start  log10 worst: start -> after GRAPE (seeds)        mean infid.   worst at 1/4, 1/2, 1 of the budget"]
    for task in TASKS:
      for m in MODES:
        for s in STARTS:
            rs = sorted((r for r in results.values() if r["task"] == task and r["mode"] == m and r["start"] == s),
                        key=lambda r: r["seed"])
            if not rs:
                continue
            lw0 = [np.log10(r["start_truth"]["worst"]) for r in rs]
            lw = [np.log10(r["truth"]["worst"]) for r in rs]
            mean = np.mean([r["truth"]["mean"] for r in rs])
            any_ = [np.mean([np.log10(r["anytime"][k]) for r in rs]) for k in sorted(rs[0]["anytime"], key=float)]
            lines.append(f"{task:8s} {m:7s} {s:5s}  {np.mean(lw0):6.2f} -> {np.mean(lw):6.2f}  ({', '.join(f'{v:.2f}' for v in lw)})"
                         f"            {mean:9.2e}   " + " / ".join(f"{v:.2f}" for v in any_))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--starts", help="folder with the BARQ start pulses UR{angle}_s{seed}.npz")
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--budget", type=float, default=300.0, help="emulated wall-clock seconds per run")
    ap.add_argument("--smoke", action="store_true")
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
        todo = [j for j in todo if j["task"] == "UR180n" and j["seed"] == 0 and j["start"] in ("rect", "barq")]
    log = open(os.path.join(a.out, "log.txt"), "a")
    print(f"{time.strftime('%H:%M:%S')} {len(todo)} runs, budget {budget:g} s emulated, {a.workers} workers",
          file=log, flush=True)
    with mp.get_context("spawn").Pool(a.workers) as pool, open(path, "a") as f:
        for r in pool.imap_unordered(_run, [(j, a.starts) for j in todo]):
            f.write(json.dumps(r) + "\n"); f.flush()
            msg = (f"{r['key']}: start {r['start_truth']['worst']:.2e} -> worst {r['truth']['worst']:.2e}, "
                   f"{r['real_s'] / 60:.1f} min" if r["ok"] else f"{r['key']}: ERROR {r['error']}")
            print(f"{time.strftime('%H:%M:%S')} {msg}", file=log, flush=True)
    text = report(load(path))
    with open(os.path.join(a.out, "report.txt"), "w") as f:
        f.write(text + "\n")
    print(f"{time.strftime('%H:%M:%S')} done\n{text}", file=log, flush=True)


if __name__ == "__main__":
    main()
