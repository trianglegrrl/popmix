#!/usr/bin/env python3
"""Render benchmark results (benchmarks/run.py JSON) as a Markdown summary.

Usage: python benchmarks/report.py benchmarks/results/<file>.json > summary.md
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def _fmt(seconds: float) -> str:
    return f"{seconds:.2f} s" if seconds < 120 else f"{seconds / 60:.1f} min"


def _reused(results: dict) -> str:
    dates = sorted({c["reference_runs"]["reused_from"] for c in results["cases"].values()
                    if "reused_from" in c["reference_runs"]})
    if not dates:
        return ""
    return (f" Reference times are from the run of {', '.join(dates)} (the reference tools are unchanged; "
            "agreement is checked against that run's outputs); popmix was run again.")


def render(results: dict) -> str:
    m, v = results["machine"], results["versions"]
    lines = [f"Measured {results['date']} on {m.get('cpu', '?')} ({m['cpus']} threads, "
             f"{m.get('memory_gb', '?')} GB); popmix {v['popmix']} ({v['popmix_commit']}), admixtools2 "
             f"{v['admixtools2']}; median of {results['repeats']} runs after a warm-up run.{_reused(results)}", "",
             "| Analysis | Reference | Reference time | popmix time (threads) | Speed-up | Peak memory, ref / popmix "
             "| Agreement | Same output at every thread count |",
             "|---|---|---|---|---|---|---|---|"]
    for name, c in results["cases"].items():
        ref = c["reference_runs"]
        threads, best = min(c["popmix"].items(), key=lambda kv: kv[1]["wall_s"]["median"])
        agree = best["agreement"]
        detail = ", ".join(f"{k} {val:.1e}" if isinstance(val, float) else f"{k} {val}"
                           for k, val in agree.items() if k != "pass")
        lines.append(
            f"| {name}: {c['analysis']} | {c['reference']} | {_fmt(ref['wall_s']['median'])} "
            f"(CV {ref['wall_s']['cv']:.1%}) | {_fmt(best['wall_s']['median'])} ({threads}; CV "
            f"{best['wall_s']['cv']:.1%}) | {ref['wall_s']['median'] / best['wall_s']['median']:.1f}x | "
            f"{ref['peak_rss_mb']['median']:.0f} / {best['peak_rss_mb']['median']:.0f} MB | "
            f"{'pass' if agree['pass'] else 'FAIL'} ({detail}) | "
            f"{'yes' if c['identical_across_threads'] else 'NO'} |")
    lines += ["", "popmix by thread count (median wall time):", ""]
    for name, c in results["cases"].items():
        per = ", ".join(f"{t}: {_fmt(r['wall_s']['median'])}" for t, r in c["popmix"].items())
        lines.append(f"- {name}: {per}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    print(render(json.loads(Path(sys.argv[1]).read_text())), end="")
