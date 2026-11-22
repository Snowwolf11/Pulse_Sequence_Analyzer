"""Quick-settings study: which GRAPE settings give the best broadband UR gates in a fixed time?

Task of the future "quick" mode: the user gives the target gate, the duration and the bandwidth; the pulse
must be smooth enough for typical NMR hardware (soft limits: amplitude rate <= 20 %/µs, phase rate <= 50 kHz,
switched on in a last stage: limits from the start blocked the optimization in pilot runs).
Here: bandwidth ±10 kHz, B1 = 10 kHz, no B1 scaling, sub-pulses of 0.5 µs, rectangular start pulse.

Scoring. Every result is judged on a fine test grid (401 offsets over the band), not on the optimization
points: main score log10(worst infidelity), also mean infidelity and the fraction of the band below 1e-3 /
1e-4 / 1e-5. A pulse that breaks the smoothness limits (> 5 % over) does not count.

Time. The quick mode will run on all worker cores (P = 11) for a wall-clock budget (20 min). Each study run
emulates that on ONE core (BLAS single-threaded): the Search tasks of a round run round-robin in slices of
a few evaluations (resumable optimize), and the emulated wall clock advances per cycle by
max(sum of slice CPU / P, longest slice), i.e. P cores, idle ones when fewer tasks than cores. The best pulse
so far on the test grid ("incumbent") is recorded at checkpoints of emulated wall time (anytime curves).

Stages (each carries the best forward; level A: all configs, short budget, one seed; later levels: the
best, longer budgets, new seeds; ranking = sum over gates of the rank of the mean log10(worst infidelity)):
    1  robustness formulation: broadband grid (density 1, 2, 4 points per 1/T) or local Taylor (density
       0.5, 1; orders 1 or 1+2), each without / with Lawson reweighting
    2  controls: full / spline (3 or 6 µs spacing), smoothness 0 / 1e-6, coarse-to-fine
    3  exploration: anchors + random combinations of starts, start pool, spread, noise correlation,
       screening, kicks (size, temperature), crossover, evaluations per candidate
    F  final: the two best recipes at the full budget on the study gates and two hold-out gates
       (other durations / axes), plus a naive baseline
A global deadline (default 14.5 h) trims the later levels and the final budget if needed.

Run:  python studies/quick_settings_study.py --out DIR [--workers 10] [--hours 14.5]
Interrupted runs resume (results.jsonl is appended; finished jobs are skipped).
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import math
import sys
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from psa.design.grape import GrapeConfig, Search, band_quality, optimize, stage_plan  # noqa: E402
from psa.design.optimize import SolverSettings  # noqa: E402
from psa.design.pulse import PulseSequence  # noqa: E402

DT, UMAX = 0.5e-6, 1e4
BAND = (-1e4, 1e4)
P = 11                                    # cores of the emulated machine (GUI: cpu_count - 1)
LIMITS = dict(max_amp_rate=20.0, max_phase_rate_khz=50.0)
ACCEPT = 1.05                             # a result may exceed the limits by 5 %
CHUNK = 20                                # evaluations per slice
POLISH = 0.2                              # last fraction of the time: the best free candidates get the limits
K_POLISH = 3                              # number of candidates polished
CHECKPOINTS = (15, 30, 60, 90, 120, 180, 240, 360, 480, 720, 900, 1200)
GATES = {
    "URx180_450": dict(axis=[1, 0, 0], angle=180.0, us=450),
    "URy90_430": dict(axis=[0, 1, 0], angle=90.0, us=430),
    "URz0_410": dict(axis=[0, 0, 1], angle=0.0, us=410),
}
HOLDOUT = {
    "URxy120_300": dict(axis=[1, 1, 0], angle=120.0, us=300),
    "URy180_600": dict(axis=[0, 1, 0], angle=180.0, us=600),
}
ALL_GATES = dict(GATES, **HOLDOUT)
CPU_PER_S = 8.0                           # CPU seconds per emulated wall second: first guess, then measured

BASE_CTRL = dict(controls="full", spline_us=0.0, smooth=0.0, coarse_us=[])
BASE_EXPL = dict(n_starts=11, start_pool="perturb", spread=0.5, corr=0.05, screen_rounds=2, screen_nfev=30,
                 keep=0.5, n_kicks=1000, kick_size=0.3, kick_temp=0.0, crossover=0.0, max_nfev=3000)
NAIVE = (dict(kind="bb", density=1, orders=[], lawson=0), BASE_CTRL,
         dict(BASE_EXPL, n_starts=1, screen_rounds=0, n_kicks=0, max_nfev=100000))


# ------------------------------------------------------------------ configs ---

def make_config(gate, form, ctrl, expl, seed):
    g = ALL_GATES[gate]
    T = g["us"] * 1e-6
    n = int(round(form["density"] * (BAND[1] - BAND[0]) * T)) + 1
    d = dict(gate="UR", target_axis=g["axis"], target_angle_deg=g["angle"], seed=int(seed),
             reweight_rounds=int(form["lawson"]), limits_late=True, **LIMITS)
    if form["kind"] == "bb":
        d.update(mode="broadband", bb_offset_min_hz=BAND[0], bb_offset_max_hz=BAND[1], bb_offset_n=n)
    else:
        d.update(mode="local", offsets_hz=[float(o) for o in np.linspace(BAND[0], BAND[1], n)], orders=list(form["orders"]))
    if ctrl["controls"] == "spline":
        d.update(controls="spline", n_spline=max(4, int(round(g["us"] / ctrl["spline_us"]))))
    d.update(smooth_weight=float(ctrl["smooth"]),
             coarse_levels=[max(4, int(round(g["us"] / c))) for c in ctrl.get("coarse_us", [])])
    d.update(n_starts=int(expl["n_starts"]), start_pool=expl["start_pool"], multistart_noise=float(expl["spread"]),
             noise_correlation_us=float(expl["corr"]) * g["us"], screen_rounds=int(expl["screen_rounds"]),
             screen_nfev=int(expl["screen_nfev"]), keep_fraction=float(expl["keep"]), n_kicks=int(expl["n_kicks"]),
             kick_size=float(expl["kick_size"]), kick_temperature=float(expl["kick_temp"]),
             crossover_rate=float(expl["crossover"]))
    cfg = GrapeConfig(**d)
    cfg.check()
    return cfg, SolverSettings(max_nfev=int(expl["max_nfev"]))


def start_pulse(gate):
    n = int(round(ALL_GATES[gate]["us"] * 1e-6 / DT))
    return PulseSequence(np.full(n, 100.0), np.zeros(n), DT, UMAX)


def rates(amp, phase):
    a = np.asarray(amp) / 100.0
    amp_rate = float(np.max(np.abs(np.diff(a)) * 100 / (DT * 1e6), initial=0.0))
    dphi = np.angle(np.exp(1j * np.radians(np.diff(phase))))
    f = np.abs(dphi) / (2 * math.pi * DT) / 1e3
    f = np.where((a[:-1] > 0.01) & (a[1:] > 0.01), f, 0.0)
    return amp_rate, float(np.max(f, initial=0.0))


def acceptable(amp, phase):
    ar, pr = rates(amp, phase)
    return ar <= LIMITS["max_amp_rate"] * ACCEPT and pr <= LIMITS["max_phase_rate_khz"] * ACCEPT


# --------------------------------------------------------------- one run ---

def run_job(job):
    """One emulated quick-mode run (single core). Returns the result record."""
    t_start, cpu_start = time.time(), time.process_time()
    rec = {k: job[k] for k in ("key", "stage", "level", "cfg_id", "gate", "seed", "budget", "spec")}
    try:
        form, ctrl, expl = job["spec"]["form"], job["spec"]["ctrl"], job["spec"]["expl"]
        cfg, st = make_config(job["gate"], form, ctrl, expl, job["seed"])
        ps0 = start_pulse(job["gate"])
        search = Search(ps0, cfg, st, parallel=P)
        budget = float(job["budget"])
        inc = {"score": math.inf, "q": None, "amp": None, "phase": None, "acceptable": False}
        best_any = {"score": math.inf, "q": None}
        wall = 0.0
        checkpoints = []
        cps = [c for c in CHECKPOINTS if c <= budget] + ([budget] if budget not in CHECKPOINTS else [])

        pool = {}                                                     # cid -> (score, w) of the best free pulses

        def consider(pulse, cid=None, w=None):
            q = band_quality(pulse, cfg, BAND[0], BAND[1])
            score = math.log10(max(q["worst infidelity"], 1e-16))
            if score < best_any["score"]:
                best_any.update(score=score, q=q)
            if cid is not None and w is not None and score < pool.get(cid, (math.inf,))[0]:
                pool[cid] = (score, w)
            if acceptable(pulse.amplitude, pulse.phase) and score < inc["score"]:
                inc.update(score=score, q=q, amp=pulse.amplitude.copy(), phase=pulse.phase.copy(), acceptable=True)

        consider(ps0)                                                 # the start pulse is the baseline

        def record(upto):
            while cps and cps[0] <= upto:
                checkpoints.append([cps.pop(0), inc["q"], best_any["q"]])

        search_end = budget * (1 - POLISH)
        while not search.done and wall < search_end:
            tasks = search.tasks()
            if not tasks:
                break
            active = [{"task": t, "state": t["state"], "left": t["budget"], "info": None} for t in tasks]
            while active and wall < search_end:
                cyc = []
                for a in list(active):
                    t = a["task"]
                    n = CHUNK if a["left"] is None else max(1, min(CHUNK, a["left"]))
                    c0 = time.process_time()
                    out, info = optimize(ps0, cfg, st, start_index=t["start_index"], state=a["state"], budget=n,
                                         start_summary=False)
                    cyc.append(time.process_time() - c0)
                    used = info["nfev"] - (a["state"]["nfev"] if a["state"] and "nfev" in a["state"] else 0)
                    a["state"], a["info"] = info["state"], info
                    if a["left"] is not None:
                        a["left"] -= max(used, 1)
                    consider(out, t["cid"], info["w"])                # not counted: cheap (no gradients)
                    if info["finished"] or (a["left"] is not None and a["left"] <= 0):
                        active.remove(a)
                        info["history"] = []
                        search.submit({"tid": t["tid"], "cid": t["cid"], "amp": out.amplitude.tolist(),
                                       "phase": out.phase.tolist(), "info": info})
                wall += max(sum(cyc) / P, max(cyc))
                record(wall)
        # polish: the best free candidates with the rate limits (the "limits" stage), in parallel
        kinds = [k for *_, k in stage_plan(cfg, ps0.N)]
        if "limits" in kinds and pool:
            picks = sorted(pool.values(), key=lambda v: v[0])[:K_POLISH]
            active = [{"state": {"stage": kinds.index("limits"), "w": w, "stage_nfev": 0, "nfev": 0}} for _, w in picks]
            while active and wall < budget:
                cyc = []
                for a in list(active):
                    c0 = time.process_time()
                    out, info = optimize(ps0, cfg, st, state=a["state"], budget=CHUNK, start_summary=False)
                    cyc.append(time.process_time() - c0)
                    a["state"] = info["state"]
                    consider(out)
                    if info["finished"]:
                        active.remove(a)
                wall += max(sum(cyc) / P, max(cyc))
                record(wall)
        record(budget)
        rec.update(ok=True, wall_emulated=wall, cpu=time.process_time() - cpu_start, real=time.time() - t_start,
                   checkpoints=checkpoints, final=inc["q"], final_acceptable=inc["acceptable"],
                   best_any=best_any["q"], score=inc["score"] if inc["acceptable"] else 0.0,
                   rates=rates(inc["amp"], inc["phase"]) if inc["amp"] is not None else None,
                   n_candidates=len(search.cands), search_log=list(search.log), search_done=search.done)
        if inc["amp"] is not None:
            np.savez_compressed(os.path.join(job["out"], "pulses", job["key"].replace("|", "_") + ".npz"),
                                amplitude=inc["amp"], phase=inc["phase"], dt=DT, umax=UMAX)
    except Exception as e:
        rec.update(ok=False, error=f"{type(e).__name__}: {e}", trace=traceback.format_exc(), score=0.0)
    return rec


# ------------------------------------------------------------ the study ---

class Study:
    def __init__(self, out, workers, hours):
        self.out, self.workers = out, workers
        os.makedirs(os.path.join(out, "pulses"), exist_ok=True)
        self.results_path = os.path.join(out, "results.jsonl")
        self.results = {}
        if os.path.exists(self.results_path):
            with open(self.results_path) as f:
                for line in f:
                    if line.strip():
                        r = json.loads(line)
                        self.results[r["key"]] = r
        meta_path = os.path.join(out, "meta.json")
        if os.path.exists(meta_path):
            self.meta = json.load(open(meta_path))
        else:
            self.meta = {"started": time.time(), "deadline": time.time() + hours * 3600}
            json.dump(self.meta, open(meta_path, "w"))
        self.deadline = self.meta["deadline"]
        self.pool = None

    def log(self, text):
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {text}"
        print(line, flush=True)
        with open(os.path.join(self.out, "progress.log"), "a") as f:
            f.write(line + "\n")

    def remaining(self):
        return self.deadline - time.time()

    # jobs
    def jobs(self, stage, level, configs, gates, seeds, budget):
        out = []
        for cid, spec in configs:
            for g in gates:
                for sd in seeds:
                    key = f"{stage}|{level}|{cid}|{g}|{sd}|{int(budget)}"
                    out.append(dict(key=key, stage=stage, level=level, cfg_id=cid, gate=g, seed=sd, budget=budget,
                                    spec=spec, out=self.out))
        return out

    def cpu_per_s(self):
        """Measured CPU seconds per emulated wall second (median of the finished runs; first guess before)."""
        vals = [r["cpu"] / r["budget"] for r in self.results.values() if r.get("ok") and r.get("budget")]
        return float(np.median(vals)) * 1.15 if len(vals) >= 5 else CPU_PER_S

    def est_cpu(self, job):
        return job["budget"] * self.cpu_per_s()

    def final_reserve(self):
        return 1200 * self.cpu_per_s() * 1.15 + 1800

    def fit(self, configs, make_jobs, reserve_s):
        """Drop the last configs until the level fits before the deadline (minus the reserve)."""
        while configs:
            jobs = make_jobs(configs)
            todo = [j for j in jobs if j["key"] not in self.results]
            waves = math.ceil(len(todo) / self.workers) if todo else 0
            longest = max((self.est_cpu(j) for j in todo), default=0)
            est = max(sum(self.est_cpu(j) for j in todo) / self.workers, longest if waves else 0) * 1.1
            if est <= self.remaining() - reserve_s or len(configs) == 1:
                return configs, jobs, est
            configs = configs[:-1]
        return configs, [], 0

    def run(self, jobs, label):
        todo = [j for j in jobs if j["key"] not in self.results]
        self.log(f"{label}: {len(jobs)} runs ({len(jobs) - len(todo)} done before), "
                 f"{self.remaining() / 3600:.1f} h left")
        if not todo:
            return
        todo.sort(key=lambda j: -j["budget"])                    # long ones first: better packing
        done = 0
        with open(self.results_path, "a") as f:
            for rec in self.pool.imap_unordered(run_job, todo):
                self.results[rec["key"]] = rec
                f.write(json.dumps(rec) + "\n"); f.flush()
                done += 1
                q = rec.get("final") or {}
                self.log(f"  {done}/{len(todo)} {rec['key']}: "
                         + (f"worst {q.get('worst infidelity', float('nan')):.2e} mean {q.get('mean infidelity', float('nan')):.2e}"
                            if rec.get("ok") and q else rec.get("error", "no acceptable pulse")))

    def rank(self, stage, level, cfg_ids, gates):
        """Config ids sorted by the sum over gates of the rank of the mean score (log10 worst infidelity)."""
        table = {}
        for cid in cfg_ids:
            for g in gates:
                vals = [r["score"] for r in self.results.values()
                        if r["stage"] == stage and r["level"] == level and r["cfg_id"] == cid and r["gate"] == g]
                table[(cid, g)] = float(np.mean(vals)) if vals else 0.0
        total = {cid: 0.0 for cid in cfg_ids}
        for g in gates:
            order = sorted(cfg_ids, key=lambda c: table[(c, g)])
            for i, c in enumerate(order):
                total[c] += i
        mean_score = {c: np.mean([table[(c, g)] for g in gates]) for c in cfg_ids}
        ranked = sorted(cfg_ids, key=lambda c: (total[c], mean_score[c]))
        summary = [{"cfg_id": c, "rank_sum": total[c], "mean_log10_worst": mean_score[c],
                    "per_gate": {g: table[(c, g)] for g in gates}} for c in ranked]
        with open(os.path.join(self.out, f"ranking_{stage}_{level}.json"), "w") as f:
            json.dump(summary, f, indent=1)
        self.log(f"ranking {stage}{level}: " + ", ".join(f"{s['cfg_id']} ({s['mean_log10_worst']:.2f})" for s in summary[:6]))
        return ranked

    def level(self, stage, level, configs, gates, seeds, budget, reserve_s):
        reserve_s = reserve_s + self.final_reserve()
        configs, jobs, est = self.fit(configs, lambda cs: self.jobs(stage, level, cs, gates, seeds, budget), reserve_s)
        self.log(f"stage {stage}{level}: {len(configs)} configs × {len(gates)} gates × {len(seeds)} seeds, "
                 f"budget {budget:.0f} s (emulated wall on {P} cores), estimated {est / 3600:.1f} h")
        self.run(jobs, f"stage {stage}{level}")
        return self.rank(stage, level, [c for c, _ in configs], gates), dict(configs)

    # stages
    def main(self):
        import multiprocessing as mp
        self.pool = mp.get_context("spawn").Pool(self.workers)
        gates = list(GATES)
        final_reserve = 0.0                   # added in level() from the measured CPU use
        try:
            # ---- stage 1: robustness formulation
            forms = []
            for c in (1, 2, 4):
                for lw in (0, 5):
                    forms.append((f"bb{c}" + ("L" if lw else ""), dict(kind="bb", density=c, orders=[], lawson=lw)))
            for c in (0.5, 1):
                for orders in ([1], [1, 2]):
                    for lw in (0, 5):
                        forms.append((f"loc{c}o{''.join(map(str, orders))}" + ("L" if lw else ""),
                                      dict(kind="local", density=c, orders=orders, lawson=lw)))
            s1 = [(fid, {"form": f, "ctrl": BASE_CTRL, "expl": BASE_EXPL}) for fid, f in forms]
            rank, specs = self.level("1", "A", s1, gates, [0], 90, final_reserve + 6 * 3600)
            top = [(c, specs[c]) for c in rank[:5]]
            rank1, specs1 = self.level("1", "B", top, gates, [1, 2], 180, final_reserve + 5 * 3600)
            best_forms = [specs1[c]["form"] for c in rank1[:2]]
            best_form_ids = rank1[:2]

            # ---- stage 2: controls and smoothness
            ctrls = []
            for cname, base in (("full", dict(controls="full", spline_us=0.0)), ("spl3", dict(controls="spline", spline_us=3.0)),
                                ("spl6", dict(controls="spline", spline_us=6.0))):
                for sm in (0.0, 1e-6):           # 1e-5 and stronger block the way here (pilot runs)
                    ctrls.append((f"{cname}_s{sm:g}", dict(base, smooth=sm, coarse_us=[])))
            for sm in (0.0,):
                ctrls.append((f"full_c2f_s{sm:g}", dict(controls="full", spline_us=0.0, smooth=sm, coarse_us=[24.0, 12.0])))
            s2 = [(f"{fid}+{cid}", {"form": f, "ctrl": c, "expl": BASE_EXPL})
                  for fid, f in zip(best_form_ids, best_forms) for cid, c in ctrls]
            rank, specs = self.level("2", "A", s2, gates, [0], 90, final_reserve + 4 * 3600)
            top = [(c, specs[c]) for c in rank[:5]]
            rank2, specs2 = self.level("2", "B", top, gates, [1, 2], 180, final_reserve + 3 * 3600)
            base_spec = specs2[rank2[0]]

            # ---- stage 3: exploration
            s3 = self.exploration_configs(base_spec, rank2[0])
            rank, specs = self.level("3", "A", s3, gates, [0], 90, final_reserve + 2.2 * 3600)
            top = [(c, specs[c]) for c in rank[:8]]
            rank, specs = self.level("3", "B", top, gates, [1], 240, final_reserve + 1.3 * 3600)
            top = [(c, specs[c]) for c in rank[:3]]
            rank3, specs3 = self.level("3", "C", top, gates, [2, 3], 480, final_reserve)

            # ---- final: full budget, study + hold-out gates, naive baseline
            budget = min(1200.0, max(300.0, (self.remaining() - 1800) / (self.cpu_per_s() * 1.15)))
            jobs = (self.jobs("F", "A", [(rank3[0], specs3[rank3[0]])], list(ALL_GATES), [4], budget)
                    + self.jobs("F", "A", [(rank3[1], specs3[rank3[1]])] if len(rank3) > 1 else [], gates, [4], budget)
                    + self.jobs("F", "A", [("naive", {"form": NAIVE[0], "ctrl": NAIVE[1], "expl": NAIVE[2]})],
                                gates[:2], [4], budget))
            self.log(f"final: budget {budget:.0f} s emulated wall (quick mode target 1200 s)")
            self.run(jobs, "final")
            json.dump({"winner": rank3[0], "spec": specs3[rank3[0]], "runner_up": rank3[1] if len(rank3) > 1 else None,
                       "runner_up_spec": specs3[rank3[1]] if len(rank3) > 1 else None, "final_budget_s": budget},
                      open(os.path.join(self.out, "winner.json"), "w"), indent=1)
            self.log("study finished")
        finally:
            self.pool.close()
            self.pool.join()

    @staticmethod
    def exploration_configs(base_spec, base_id, n_random=21, seed=12345):
        rng = np.random.default_rng(seed)
        anchors = {
            "anchor_base": dict(BASE_EXPL),
            "anchor_explore": dict(n_starts=24, start_pool="mix", spread=1.0, corr=0.1, screen_rounds=3, screen_nfev=25,
                                   keep=0.5, n_kicks=1000, kick_size=0.5, kick_temp=0.3, crossover=0.0, max_nfev=3000),
            "anchor_kicks_only": dict(BASE_EXPL, n_starts=1, screen_rounds=0, kick_size=0.3, kick_temp=0.1),
        }
        space = dict(n_starts=[1, 4, 11, 22, 33], start_pool=["perturb", "mix", "random"], spread=[0.2, 0.5, 1.0],
                     corr=[0.0, 0.02, 0.05, 0.1], screen_rounds=[0, 1, 2, 3], screen_nfev=[15, 30, 60], keep=[0.33, 0.5],
                     n_kicks=[0, 1000, 1000], kick_size=[0.1, 0.3, 0.6], kick_temp=[0.0, 0.1, 0.3],
                     crossover=[0.0, 0.0, 0.3], max_nfev=[1000, 3000, 10000])
        # balanced random design: every value of every factor appears about equally often
        cols = {}
        for k, vals in space.items():
            col = [vals[i % len(vals)] for i in range(n_random)]
            rng.shuffle(col)
            cols[k] = col
        out = [(f"{base_id}|{name}", dict(base_spec, expl=e)) for name, e in anchors.items()]
        for i in range(n_random):
            e = {k: cols[k][i] for k in space}
            e = {k: (v.item() if hasattr(v, "item") else v) for k, v in e.items()}
            out.append((f"{base_id}|x{i:02d}", dict(base_spec, expl=e)))
        return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--hours", type=float, default=14.5)
    ap.add_argument("--smoke", action="store_true", help="tiny budgets, to check the pipeline")
    args = ap.parse_args()
    if args.smoke:
        CPU_PER_S = 1.0
        Study.final_reserve = lambda self: 0.0
        _orig = Study.jobs

        def _small(self, stage, level, configs, gates, seeds, budget):
            return _orig(self, stage, level, configs[:2], gates[:1], seeds[:1], 15.0)
        Study.jobs = _small
        Study.exploration_configs = staticmethod(lambda b, i, **k: [(f"{i}|anchor_base", dict(b, expl=BASE_EXPL)),
                                                                     (f"{i}|anchor_kicks", dict(b, expl=dict(BASE_EXPL, n_starts=1)))])
    try:
        os.nice(10)                       # leave the machine responsive for other work
    except OSError:
        pass
    Study(args.out, args.workers, args.hours).main()
