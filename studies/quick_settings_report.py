"""Summary of a quick-settings study (studies/quick_settings_study.py): rankings per stage, factor effects in
the exploration stage, final results, anytime curves and fidelity-vs-offset plots of the final pulses.

Run:  python studies/quick_settings_report.py DIR      (writes DIR/report/summary.md and PNG plots;
works on a running study, too)
"""
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)


def load(out):
    res = []
    with open(os.path.join(out, "results.jsonl")) as f:
        for line in f:
            if line.strip():
                res.append(json.loads(line))
    return res


def fmt(x):
    return "–" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.2e}"


def level_table(res, stage, level):
    rows = defaultdict(lambda: defaultdict(list))
    for r in res:
        if r["stage"] == stage and r["level"] == level:
            rows[r["cfg_id"]][r["gate"]].append(r)
    if not rows:
        return None
    gates = sorted({g for v in rows.values() for g in v})
    table = []
    for cid, per in rows.items():
        scores = {g: np.mean([x["score"] for x in per[g]]) for g in gates if per[g]}
        worst = {g: np.median([(x.get("final") or {}).get("worst infidelity", 1.0) for x in per[g]]) for g in gates if per[g]}
        mean = {g: np.median([(x.get("final") or {}).get("mean infidelity", 1.0) for x in per[g]]) for g in gates if per[g]}
        table.append((float(np.mean(list(scores.values()))), cid, worst, mean, sum(len(v) for v in per.values())))
    table.sort()
    return gates, table


def factor_effects(res, stage="3", level="A"):
    """Mean score (log10 worst infidelity) per value of every exploration factor."""
    eff = defaultdict(lambda: defaultdict(list))
    for r in res:
        if r["stage"] == stage and r["level"] == level and r.get("ok"):
            for k, v in r["spec"]["expl"].items():
                eff[k][json.dumps(v)].append(r["score"])
    return {k: {v: (float(np.mean(s)), len(s)) for v, s in sorted(d.items())} for k, d in eff.items()}


def plots(out, res, report):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    from psa.design.grape import GrapeConfig, band_infidelity
    from psa.design.pulse import PulseSequence
    import quick_settings_study as Q
    files = []
    final = [r for r in res if r["stage"] == "F" and r.get("ok")]
    if not final:                                     # study still running: use the latest level
        last = max(res, key=lambda r: (r["stage"], r["level"]), default=None)
        final = [r for r in res if last and r["stage"] == last["stage"] and r["level"] == last["level"] and r.get("ok")]
    if not final:
        return files
    gates = sorted({r["gate"] for r in final})
    # anytime curves
    fig, axes = plt.subplots(1, len(gates), figsize=(4.2 * len(gates), 3.4), squeeze=False)
    for ax, g in zip(axes[0], gates):
        for r in sorted([x for x in final if x["gate"] == g], key=lambda x: x["cfg_id"]):
            free = [(c[0], c[2]["worst infidelity"]) for c in r["checkpoints"] if len(c) > 2 and c[2]]
            acc = [(c[0], c[1]["worst infidelity"]) for c in r["checkpoints"] if c[1]]
            line = None
            if free:
                t, w = zip(*free)
                line, = ax.semilogy(t, w, ls="--", lw=1, label=r["cfg_id"].split("|")[-1] + " (free)")
            if acc:
                t, w = zip(*acc)
                ax.semilogy(t, w, marker="o", ms=3, color=line.get_color() if line else None,
                            label=r["cfg_id"].split("|")[-1] + " (acceptable)")
        ax.set_title(g); ax.set_xlabel("emulated wall time on 11 cores [s]"); ax.grid(alpha=0.3)
    axes[0][0].set_ylabel("worst infidelity over ±10 kHz")
    axes[0][-1].legend(fontsize=7)
    fig.tight_layout()
    path = os.path.join(report, "anytime.png"); fig.savefig(path, dpi=130); plt.close(fig); files.append(path)
    # infidelity vs offset of the final pulses
    fig, axes = plt.subplots(1, len(gates), figsize=(4.2 * len(gates), 3.4), squeeze=False)
    offs = np.linspace(-12e3, 12e3, 601)
    for ax, g in zip(axes[0], gates):
        spec = Q.ALL_GATES[g]
        cfg = GrapeConfig(mode="broadband", gate="UR", target_axis=spec["axis"], target_angle_deg=spec["angle"])
        for r in sorted([x for x in final if x["gate"] == g], key=lambda x: x["cfg_id"]):
            path = os.path.join(out, "pulses", r["key"].replace("|", "_") + ".npz")
            if not os.path.exists(path):
                continue
            d = np.load(path)
            ps = PulseSequence(d["amplitude"], d["phase"], float(d["dt"]), float(d["umax"]))
            inf = np.maximum(band_infidelity(ps, cfg, offs), 1e-16)
            ax.semilogy(offs / 1e3, inf, lw=1, label=r["cfg_id"].split("|")[-1])
        ax.axvspan(-10, 10, color="0.9", zorder=0)
        ax.set_title(g); ax.set_xlabel("offset [kHz]"); ax.grid(alpha=0.3)
    axes[0][0].set_ylabel("infidelity")
    axes[0][-1].legend(fontsize=7)
    fig.tight_layout()
    path = os.path.join(report, "offset_profiles.png"); fig.savefig(path, dpi=130); plt.close(fig); files.append(path)
    return files


def main(out):
    res = load(out)
    report = os.path.join(out, "report")
    os.makedirs(report, exist_ok=True)
    lines = [f"# Quick-settings study – {len(res)} runs", ""]
    if os.path.exists(os.path.join(out, "winner.json")):
        w = json.load(open(os.path.join(out, "winner.json")))
        lines += ["## Winner", "", f"`{w['winner']}` (final budget {w['final_budget_s']:.0f} s)", "",
                  "```json", json.dumps(w["spec"], indent=1), "```", ""]
    for stage, level in [("1", "A"), ("1", "B"), ("2", "A"), ("2", "B"), ("3", "A"), ("3", "B"), ("3", "C"), ("F", "A")]:
        t = level_table(res, stage, level)
        if t is None:
            continue
        gates, table = t
        lines += [f"## Stage {stage}{level}", "", "| config | mean log10 worst | " + " | ".join(f"{g} worst / mean" for g in gates) + " | runs |",
                  "|---|---|" + "---|" * len(gates) + "---|"]
        for score, cid, worst, mean, n in table:
            lines.append(f"| {cid} | {score:.2f} | " + " | ".join(f"{fmt(worst.get(g))} / {fmt(mean.get(g))}" for g in gates) + f" | {n} |")
        lines.append("")
    eff = factor_effects(res)
    if eff:
        lines += ["## Exploration factors (stage 3A: mean log10 worst infidelity per value; lower is better)", ""]
        for k, d in eff.items():
            lines.append(f"- **{k}**: " + ", ".join(f"{v} → {m:.2f} ({n})" for v, (m, n) in d.items()))
        lines.append("")
    errors = [r for r in res if not r.get("ok")]
    if errors:
        lines += ["## Errors", ""] + [f"- {r['key']}: {r.get('error')}" for r in errors[:20]] + [""]
    files = plots(out, res, report)
    lines += [f"![{os.path.basename(f)}]({os.path.basename(f)})" for f in files]
    with open(os.path.join(report, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(os.path.join(report, "summary.md"))


if __name__ == "__main__":
    main(sys.argv[1])
