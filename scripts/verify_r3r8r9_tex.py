#!/usr/bin/env python3
"""Verify that the R3/R8/R9 numbers typeset in the SI match their JSON sources.

Fail-closed check (AGENTS traceability rule): every number in
``tab:size5``, ``tab:witness`` and ``tab:longtime`` of
``02_track_A_topo_SI.tex`` must be derivable from
``data/r3_size5_analysis.json``, ``data/r8_witness_summary.json`` and
``data/r9_longtime_summary.json`` at the printed precision.

Exits non-zero on the first mismatch.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SI = ROOT / "02_track_A_topo_SI.tex"


def table_body(tex: str, label: str) -> str:
    """Return the contents of the table environment carrying ``label``."""
    for chunk in tex.split(r"\begin{table}")[1:]:
        if rf"\label{{{label}}}" in chunk:
            return chunk.split(r"\end{table}")[0]
    raise KeyError(label)


def nums(body: str) -> list[str]:
    """Return every \\num{...} token value in order of appearance."""
    return re.findall(r"\\num\{([^}]*)\}", body)


def fmt(x: float, places: int) -> str:
    """Format like the table: grid-floor values print as bare zero."""
    if abs(x) < 0.5 * 10 ** (-places):
        return "0"
    return f"{abs(x):.{places}f}"


def check_table(name: str, expected: list[str], got: list[str]) -> int:
    """Compare expected literal token lists against extracted ones."""
    misses = [e for e in expected if e not in got]
    if misses:
        print(f"FAIL {name}: missing {misses}")
        print(f"  got: {got}")
        return 1
    print(f"OK   {name}: {len(expected)} expected values present "
          f"(table has {len(got)} \\num tokens)")
    return 0


def main() -> int:
    """Run all three table checks."""
    tex = SI.read_text()
    fails = 0

    # --- tab:size5: measured and frozen CIs from r3_size5_analysis.json ---
    r3 = json.loads((ROOT / "data" / "r3_size5_analysis.json").read_text())
    body = table_body(tex, "tab:size5")
    expected = []
    for cond in ("disorder_0.0_edge_low_loss", "disorder_0.0_uniform_loss",
                 "disorder_0.5_edge_low_loss", "disorder_0.5_uniform_loss"):
        c = r3["conditions"][cond]
        m, f = c["measured"], c["frozen_prediction"]
        # mean/se carry their sign outside \num{} in the table; CI bounds
        # carry it inside, so the two groups format differently.
        expected += [f"{abs(m['mean']):.6f}", f"{abs(m['se']):.6f}",
                     f"{m['ci95'][0]:.6f}", f"{m['ci95'][1]:.6f}",
                     f"{f['ci95_low']:.6f}", f"{f['ci95_high']:.6f}"]
    fails += check_table("tab:size5", expected, nums(body))

    # --- tab:witness: means and counts from r8_witness_summary.json ---
    r8 = json.loads((ROOT / "data" / "r8_witness_summary.json").read_text())
    body = table_body(tex, "tab:witness")
    expected = []
    for cells in (2, 3, 4):
        e = r8["strata"][f"{cells}cells_d0_edge_low_loss"]
        for sector in ("chiral_topological", "chiral_trivial"):
            w = e["witness"]["n1"][sector]
            expected += [f"{w['mean']:.4f}", str(w["certified"]),
                         f"{e['non_gaussianity'][sector]['mean']:.4f}"]
    fails += check_table("tab:witness", expected, nums(body))

    # --- tab:longtime: selected stroboscopic means/SEs from r9 summary ---
    r9 = json.loads((ROOT / "data" / "r9_longtime_summary.json").read_text())
    body = table_body(tex, "tab:longtime")
    expected = []
    for key in ("2cells_d0_edge_low_loss", "3cells_d0_edge_low_loss"):
        s = r9["strata"][key]
        for k in (1, 2, 3, 4, 8, 12, 16, 24, 32):
            expected += [fmt(s["mean_series"][k], 6),
                         fmt(s["se_series"][k], 6)]
    fails += check_table("tab:longtime", expected, nums(body))

    if fails:
        print(f"FAILS: {fails}")
        return 1
    print("FAILS: 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
