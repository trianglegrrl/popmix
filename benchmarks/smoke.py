"""Time popmix's main commands on the small dataset shipped in the repository.

The CI subset of the benchmarks (TASKS.md E6): no reference tool, so it
shows only how popmix's own times move from commit to commit. The input is
tests/golden/dates/family_chr20_22 (54 samples, 21,162 SNPs), so the times
are mostly startup and reading; `run.py` has the full comparison.

    python benchmarks/smoke.py [--repeats 3] [--out summary.md]
"""

from __future__ import annotations

import argparse
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "golden" / "dates" / "family_chr20_22"


def commands(out: Path) -> dict[str, list[str]]:
    """Each command's popmix arguments, writing under `out`."""
    i = str(DATA)
    return {
        "extract-f2": ["extract-f2", "-i", i, "-o", str(out / "f2"), "--populations", "French,Yoruba,Simulation.1"],
        "qpdstat": ["qpdstat", "-i", i, "--pop1", "French", "--pop2", "Yoruba", "--pop3", "Simulation.1",
                    "--pop4", "Simulation"],
        "pca": ["pca", "-i", i, "-o", str(out / "pca")],
        "kinship": ["kinship", "-i", i, "--method", "king", "-o", str(out / "king.tsv")],
        "dates": ["dates", "-i", i, "--admixed", "Simulation.1", "--sources", "French,Yoruba"],
    }


def run(repeats: int) -> dict[str, list[float]]:
    """Wall times of each command, `repeats` times after one warm-up run."""
    times: dict[str, list[float]] = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, args in commands(Path(tmp)).items():
            argv = ["popmix", "--threads", "2", *args]
            times[name] = []
            for k in range(repeats + 1):
                start = time.perf_counter()
                subprocess.run(argv, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                if k:
                    times[name].append(time.perf_counter() - start)
    return times


def render(times: dict[str, list[float]]) -> str:
    """A Markdown table of the median and range of each command's times."""
    rows = ["| Command | Median | Min | Max |", "|---|---|---|---|"]
    for name, ts in times.items():
        rows.append(f"| {name} | {statistics.median(ts):.2f} s | {min(ts):.2f} s | {max(ts):.2f} s |")
    return f"popmix on {DATA.name} (54 samples, 21,162 SNPs), 2 threads:\n\n" + "\n".join(rows) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--out", type=Path, default=None, help="Also append the table to this file")
    args = parser.parse_args()
    table = render(run(args.repeats))
    sys.stdout.write(table)
    if args.out is not None:
        with args.out.open("a") as f:
            f.write(table)


if __name__ == "__main__":
    main()
