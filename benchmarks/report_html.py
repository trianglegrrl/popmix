#!/usr/bin/env python3
"""Render benchmark results as one self-contained HTML page (TASKS.md E6).

The page (report_template.html) draws its tables and charts from the results
of run.py and, if given, scaling.py, embedded as JSON.

Usage: python benchmarks/report_html.py results/<host>-<date>.json
                                        [--scaling results/scaling-<host>-<date>.json] [-o report.html]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

TEMPLATE = Path(__file__).parent / "report_template.html"


def render(results: dict, scaling: dict | None = None) -> str:
    data = {**results, "scaling": scaling}
    # "</" cannot end the script element from inside the JSON
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    return TEMPLATE.read_text().replace("/*DATA*/null", payload, 1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", type=Path)
    parser.add_argument("--scaling", type=Path)
    parser.add_argument("-o", "--out", type=Path, default=Path("benchmark-report.html"))
    args = parser.parse_args()
    scaling = json.loads(args.scaling.read_text()) if args.scaling else None
    args.out.write_text(render(json.loads(args.results.read_text()), scaling))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
