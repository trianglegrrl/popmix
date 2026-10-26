#!/usr/bin/env python3
"""Benchmark popmix against the reference tools (TASKS.md E1-E3).

For each case: popmix at each thread count and the reference tool, each run
`--repeats` times after a warm-up run, on the same input; wall time, CPU time
and peak memory per run, with median, IQR and coefficient of variation;
agreement of the outputs at the CLAUDE.md tolerances; and whether popmix's
output is byte-identical at every thread count.

Usage: python benchmarks/run.py [--cases extract_f2,qpadm] [--repeats 5]
                                [--threads 1,4,16] [--out results.json]
                                [--reuse-reference earlier.json]
Results go to benchmarks/results/<host>-<date>.json unless --out is given.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import shutil
import socket
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from cases import CACHE, CASES  # noqa: E402
from harness import machine, measure  # noqa: E402


def _digest(directory: Path) -> str:
    h = hashlib.sha256()
    for path in sorted(p for p in directory.rglob("*") if p.is_file() and p.suffix not in (".par", ".R", ".log")):
        h.update(path.name.encode())
        h.update(path.read_bytes())
    return h.hexdigest()


def _fresh(directory: Path) -> Path:
    shutil.rmtree(directory, ignore_errors=True)
    directory.mkdir(parents=True)
    return directory


def reused_reference(earlier: dict, name: str, ref_dir: Path) -> dict:
    """The reference timing of case `name` from an earlier report, whose
    outputs are still in `ref_dir`, marked with the earlier date."""
    case = earlier.get("cases", {}).get(name)
    if case is None:
        raise SystemExit(f"the earlier report has no case {name}")
    if not ref_dir.is_dir() or not any(ref_dir.iterdir()):
        raise SystemExit(f"{name}: no reference outputs in {ref_dir}")
    return {**case["reference_runs"], "reused_from": earlier["date"]}


def versions() -> dict:
    def out(cmd: list[str]) -> str:
        try:
            return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()
        except (OSError, subprocess.TimeoutExpired):
            return "unavailable"
    import popmix

    return {"popmix": popmix.__version__,
            "popmix_commit": out(["git", "-C", str(Path(__file__).parent), "rev-parse", "--short", "HEAD"]),
            "admixtools2": out(["/usr/bin/Rscript", "-e", "cat(as.character(packageVersion('admixtools')))"]),
            "R": out(["/usr/bin/Rscript", "-e", "cat(R.version.string)"])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cases", default=",".join(CASES))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--threads", default="1,4,16")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--reuse-reference", type=Path, metavar="REPORT",
                        help="take the reference timings from this earlier report and its outputs from the "
                             "run directories, and run only popmix")
    args = parser.parse_args()
    earlier = json.loads(args.reuse_reference.read_text()) if args.reuse_reference else None
    threads = [int(t) for t in args.threads.split(",")]
    today = datetime.date.today().isoformat()
    out = args.out or Path(__file__).parent / "results" / f"{socket.gethostname()}-{today}.json"
    report = {"date": today, "machine": machine(), "versions": versions(), "repeats": args.repeats, "cases": {}}

    for name in args.cases.split(","):
        case = CASES[name]
        source = case.prepare()
        work = CACHE / "runs" / name
        if earlier is not None:
            ref_dir = work / "reference"
            ref_runs = reused_reference(earlier, name, ref_dir)
        else:
            print(f"{name}: reference ({case.reference})", flush=True)
            ref_dir = _fresh(work / "reference")
            ref_runs = measure(case.ref(source, ref_dir), ref_dir, args.repeats, log=ref_dir / "run.log").summary()
        result = {"analysis": case.analysis, "reference": case.reference,
                  "reference_runs": ref_runs, "popmix": {}}
        digests = {}
        for t in threads:
            print(f"{name}: popmix, {t} threads", flush=True)
            ours = _fresh(work / f"popmix-t{t}")
            m = measure(case.popmix(source, ours, t), ours, args.repeats, log=ours / "run.log")
            result["popmix"][str(t)] = m.summary()
            result["popmix"][str(t)]["agreement"] = case.agree(ours, ref_dir)
            digests[t] = _digest(ours)
        result["identical_across_threads"] = len(set(digests.values())) == 1
        report["cases"][name] = result
        _write(out, report)
        best = min(result["popmix"].values(), key=lambda r: r["wall_s"]["median"])
        print(f"  reference {ref_runs['wall_s']['median']:.2f}s, popmix best {best['wall_s']['median']:.2f}s, "
              f"agreement {best['agreement']}, identical across threads: {result['identical_across_threads']}",
              flush=True)

    _write(out, report)
    print(f"wrote {out}")


def _plain(value):
    """numpy scalars as Python ones, for JSON."""
    return value.item() if hasattr(value, "item") else str(value)


def _write(out: Path, report: dict) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1, default=_plain) + "\n")


if __name__ == "__main__":
    main()
