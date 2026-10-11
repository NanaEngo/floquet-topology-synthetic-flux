#!/usr/bin/env python3
"""Audit and summarize the chiral Floquet production outputs.

The independent unit is the matched disorder realization.  Each contrast is
computed within (size, loss profile, disorder level, sample), before any
aggregation.  This script is descriptive and does not certify topological
protection or causal mechanism.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SECTORS = ("chiral_topological", "chiral_trivial")
PROFILES = ("edge_low_loss", "uniform_loss")
DISORDER_LEVELS = (0.0, 0.5, 1.0, 1.5)


def _dedupe_by_cells(files):
    by_cells = {}
    for path in files:
        data = json.loads(path.read_text())
        cells = int(data["design"]["n_cells"])
        n_mc = int(data["design"]["n_mc"])
        prev = by_cells.get(cells)
        if prev is None or n_mc > int(json.loads(prev.read_text())["design"]["n_mc"]):
            by_cells[cells] = path
    return [by_cells[c] for c in sorted(by_cells)]


def finite_float(value):
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def paired_rows(rows):
    groups = {}
    for row in rows:
        key = (row["loss_profile"], float(row["disorder_level"]), int(row["sample"]), int(row["draw_seed"]))
        groups.setdefault(key, {})[row["sector"]] = row
    if any(set(group) != set(SECTORS) for group in groups.values()):
        raise ValueError("incomplete sector pairing detected")
    return groups


def mean_ci(values):
    x = np.asarray(values, dtype=float)
    n = len(x)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=1)) if n > 1 else 0.0
    se = sd / math.sqrt(n) if n else float("nan")
    # Normal 95% interval keeps the summary dependency-light; the raw paired
    # values are retained in the production JSON for any later bootstrap.
    half = 1.96 * se
    return {"n": n, "mean": mean, "sd": sd, "se": se,
            "ci95_low": mean - half, "ci95_high": mean + half,
            "positive_fraction": float(np.mean(x > 0)) if n else float("nan")}


def load(path):
    data = json.loads(path.read_text())
    rows = data.get("records", [])
    if not rows:
        raise ValueError(f"no records: {path}")
    expected = 2 * len(PROFILES) * len(DISORDER_LEVELS) * int(data["design"]["n_mc"])
    if len(rows) != expected:
        raise ValueError(f"{path.name}: expected {expected} records, found {len(rows)}")
    for row in rows:
        for value in row.values():
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError(f"non-finite value in {path.name}")
    groups = paired_rows(rows)
    expected_groups = len(PROFILES) * len(DISORDER_LEVELS) * int(data["design"]["n_mc"])
    if len(groups) != expected_groups:
        raise ValueError(f"{path.name}: expected {expected_groups} paired units, found {len(groups)}")
    return data, groups


def summarize(path):
    data, groups = load(path)
    cells = int(data["design"]["n_cells"])
    out = []
    for profile in PROFILES:
        for disorder in DISORDER_LEVELS:
            values = []
            center_values = []
            edge_pop_values = []
            bulk_pop_values = []
            for (p, d, _sample, _seed), group in sorted(groups.items()):
                if p != profile or d != disorder:
                    continue
                top = group["chiral_topological"]
                trivial = group["chiral_trivial"]
                values.append(top["negativity_edge"] - trivial["negativity_edge"])
                center_values.append(top["negativity_center"] - trivial["negativity_center"])
                edge_pop_values.append(top["edge_population"] - trivial["edge_population"])
                bulk_pop_values.append(top["bulk_population"] - trivial["bulk_population"])
            s = mean_ci(values)
            s.update({"n_cells": cells, "loss_profile": profile,
                      "disorder_level": disorder,
                      "observable": "negativity_edge_topological_minus_trivial",
                      "center_negativity": mean_ci(center_values),
                      "edge_population": mean_ci(edge_pop_values),
                      "bulk_population": mean_ci(bulk_pop_values)})
            out.append(s)
    return {"input": str(path.relative_to(ROOT)), "n_cells": cells,
            "n_records": len(data["records"]), "n_paired_units": len(groups),
            "design": data["design"], "summaries": out}


def _collective_panel_data(collective_analysis: dict) -> dict:
    """Extract per-truncation (profile, disorder) paired contrasts from a
    collective-mode analysis JSON (keys ``paired_contrasts.<observable>``)."""
    # keep the analysis-JSON keys verbatim ("edge_low_loss_d0.0", ...)
    return {observable: dict(table)
            for observable, table in collective_analysis.get("paired_contrasts", {}).items()}


def _fivecell_points(path):
    """Out-of-sample five-cell contrasts from the frozen-prediction scoring.

    Returns {(profile, disorder): (mean, se)} for the two measured
    disorder levels (0.0, 0.5); keys look like "disorder_0.0_edge_low_loss".
    """
    data = json.loads(path.read_text())
    points = {}
    for key, cond in data.get("conditions", {}).items():
        parts = key.split("_")
        points[("_".join(parts[2:]), float(parts[1]))] = (
            cond["measured"]["mean"], cond["measured"]["se"])
    return points


def plot(summaries, output, collective_analyses=None, fivecell=None):
    import matplotlib.pyplot as plt
    plt.style.use(["seaborn-v0_8-ticks", "seaborn-v0_8-colorblind"])
    # file-summaries carry per-condition rows under ``summaries``; flatten once.
    if summaries and "summaries" in summaries[0]:
        summaries = [item for fs in summaries for item in fs["summaries"]]
    sizes = sorted({x["n_cells"] for x in summaries})
    n_panels = len(PROFILES) + (1 if collective_analyses else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(11.5, 4.3))
    if n_panels == 1:
        axes = [axes]
    colors = {0.0: "#1b9e77", 0.5: "#d95f02", 1.0: "#7570b3", 1.5: "#e7298a"}
    for ax, profile, letter in zip(axes[: len(PROFILES)], PROFILES, ["(a)", "(b)"]):
        ax.text(0.02, 0.98, letter, transform=ax.transAxes,
                va="top", ha="left", fontsize=11)
        for disorder in DISORDER_LEVELS:
            xs, ys, lo, hi = [], [], [], []
            for size in sizes:
                match = next((x for x in summaries if x["n_cells"] == size and x["loss_profile"] == profile and x["disorder_level"] == disorder), None)
                if match is None:
                    continue
                xs.append(size); ys.append(match["mean"])
                lo.append(match["ci95_low"]); hi.append(match["ci95_high"])
            if xs:
                ax.errorbar(xs, ys, yerr=[np.asarray(ys)-np.asarray(lo), np.asarray(hi)-np.asarray(ys)],
                            marker="o", linewidth=1.5, capsize=3, color=colors[disorder],
                            label=f"disorder={disorder:g}")
        # Out-of-sample five-cell point (frozen-prediction scoring, N_MC=300,
        # disorders 0.0/0.5 only): star markers, no connecting line.
        if fivecell:
            for disorder in (0.0, 0.5):
                key = (profile, disorder)
                if key not in fivecell:
                    continue
                mean, se = fivecell[key]
                first = (ax is axes[1] and disorder == 0.0)
                ax.errorbar([5], [mean], yerr=[1.96 * se], marker="*",
                            markersize=10, linestyle="None", capsize=3,
                            color=colors[disorder],
                            markeredgecolor=colors[disorder],
                            label="5 cells, out-of-sample ($N_{\\rm MC}=300$)" if first else "_nolegend_")
        ax.axhline(0, color="0.25", linewidth=0.8, linestyle="--")
        ax.set_xlabel("Number of unit cells")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel(r"$\Delta\mathcal{N}_{\rm edge}$ (topological $-$ trivial)")
    axes[0].set_ylim(-0.006, 0.022)
    axes[1].set_ylim(-0.002, 0.002)
    axes[min(len(PROFILES), n_panels) - 1].legend(frameon=False, fontsize=8)
    # --- third panel: local vs collective at two cells, both truncations ---
    if collective_analyses:
        ax = axes[len(PROFILES)]
        local_d0 = {}; local_d35 = {}; collective_d0 = {}; collective_d35 = {}
        for truncation, data in collective_analyses.items():
            local = data.get("local_site0", {})
            coll = data.get("collective_edge", {})
            local_d0[truncation] = local.get("edge_low_loss_d0.0")
            local_d35[truncation] = local.get("edge_low_loss_d0.35")
            collective_d0[truncation] = coll.get("edge_low_loss_d0.0")
            collective_d35[truncation] = coll.get("edge_low_loss_d0.35")
        xs = sorted(collective_analyses)
        xloc0 = [x - 0.05 for x in xs]   # d = 0 series
        xloc35 = [x + 0.05 for x in xs]  # d = 0.35 series
        for xoff, series, marker, fc, label in [
            (xloc0, local_d0, "o", "#1b9e77", "local site-0, $d=0$"),
            (xloc0, collective_d0, "s", "white", "collective edge mode, $d=0$"),
            (xloc35, local_d35, "o", "white", "local site-0, $d=0.35$"),
            (xloc35, collective_d35, "s", "white", "collective edge mode, $d=0.35$"),
        ]:
            ys = [series[x]["mean"] for x in xs]
            errs = [[series[x].get("se", 0)] for x in xs]
            color = "#1b9e77" if "local" in label else ("#d62728" if "$d=0$" in label else "#ff7f0e")
            ax.errorbar(xoff, ys, yerr=errs, marker=marker, linestyle="-",
                        linewidth=1.5, capsize=3, markersize=7,
                        color=color, markerfacecolor=fc, markeredgecolor=color,
                        label=label)
        ax.axhline(0, color="0.25", linewidth=0.8, linestyle="--")
        ax.text(0.02, 0.98, "(c)", transform=ax.transAxes,
                va="top", ha="left", fontsize=11)
        ax.set_xlabel(r"Fock truncation $N_{\max}$")
        ax.set_xticks(xs)
        ax.set_ylim(-0.003, 0.022)
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, fontsize=7, loc="upper right", ncol=1)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, default=ROOT / "data" / "chiral_production")
    p.add_argument("--output", type=Path, default=ROOT / "data" / "chiral_production_analysis.json")
    p.add_argument("--csv", type=Path, default=ROOT / "data" / "chiral_production_summary.csv")
    p.add_argument("--figure", type=Path, default=ROOT / "figures" / "chiral_production_contrast.png")
    p.add_argument("--collective-nmax2", type=Path, default=ROOT / "data" / "collective_wigner_analysis.json")
    p.add_argument("--collective-nmax3", type=Path, default=ROOT / "data" / "collective_wigner_analysis_nmax3.json")
    p.add_argument("--fivecell", type=Path, default=ROOT / "data" / "r3_size5_analysis.json")
    args = p.parse_args()
    files = sorted(args.input_dir.glob("chiral_floquet_*cells_NMC*.json"))
    # The size-5 out-of-sample run carries only two disorder levels and is
    # scored by analyze_r3_size5.py; this summary covers the four-level
    # production protocol only.
    files = [f for f in files
             if tuple(json.loads(f.read_text())["design"].get("disorder_levels", DISORDER_LEVELS)) == DISORDER_LEVELS]
    if not files:
        raise SystemExit("no production files found")
    # If multiple files share a cell size (e.g. an N_MC-superset rerun), keep
    # only the largest N_MC so the summary does not double-count a size.
    files = _dedupe_by_cells(files)
    summaries = [summarize(path) for path in files]
    output = {"status": "PRODUCTION_OUTPUTS_AUDITED", "inputs": [x["input"] for x in summaries],
              "files": summaries,
              "claim_boundary": ["descriptive matched dynamics only", "no topological protection claim", "no experimental or metrological claim"]}
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    collective_analyses = {}
    if args.collective_nmax2.is_file():
        collective_analyses[2] = _collective_panel_data(json.loads(args.collective_nmax2.read_text()))
    if args.collective_nmax3.is_file():
        collective_analyses[3] = _collective_panel_data(json.loads(args.collective_nmax3.read_text()))
    fivecell = _fivecell_points(args.fivecell) if args.fivecell.is_file() else None
    plot(summaries, args.figure, collective_analyses=collective_analyses or None, fivecell=fivecell)
    rows = []
    for file_summary in summaries:
        for item in file_summary["summaries"]:
            rows.append({k: v for k, v in item.items() if not isinstance(v, dict)})
    import csv
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"status": output["status"], "files": len(files), "summary_rows": len(rows), "figure": str(args.figure)}, indent=2))


if __name__ == "__main__":
    main()
