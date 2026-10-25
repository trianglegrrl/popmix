#!/usr/bin/env python3
"""Generate the ROLLOFF (ADMIXTOOLS 7.0.2) reference for popmix's rolloff.

Usage: python tests/golden/rolloff/generate_rolloff.py [config...]

Runs rolloff built to print full precision (rolloff-hiprec, from
build_rolloff.sh) and as released (rolloff-stock), checks they agree to the
stock build's printed precision, then follows ADMIXTOOLS' expfit.sh: fits
C + A exp(-m d) to the correlation column (4) from 0.5 to 100 cM with
rexpfit.r's method (DEoptim start, seed 1, then nls; rexpfit_full.R, checked
against rexpfit.r itself to the 8 digits it prints) for the full output and
each chromosome-left-out output, and runs dowtjack with each chromosome
weighted by its SNPs in the .snp file (wtjack.pl's weights).

Writes tests/golden/rolloff/v7.0.2/<config>/rolloff.json with the params,
the output bins (d cM, weighted Z, regression, correlation, pair count),
the jackknife bins per chromosome (pair count on that chromosome,
regression, correlation, pairs without it), the fits (C, A, m) and the
jackknife (mean, standard error).
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AADR = "/references/aadr/v62.0_1240k_public"
SMALL = str(REPO / "tests/golden/dates/family_chr20_22")
BIN = Path(os.environ.get("ROLLOFF_DIR", Path.home() / ".local/share/popmix-references/admixtools-7.0.2"))

PARAMS = {"binsize": 0.001, "maxdis": 0.2, "runmode": 1, "chithresh": 6.0, "zdipcorrmode": "YES",
          "checkmap": "YES", "seed": 1, "numchrom": 22, "jackknife": "YES"}
CONFIGS = {
    "family": {"prefix": SMALL, "admixpop": "Simulation.1", "refs": ["French", "Yoruba"], "chromosomes": None},
    "aadr_asw": {"prefix": AADR, "admixpop": "ASW.DG", "refs": ["YRI.DG", "CEU.DG"],
                     "chromosomes": ["18", "19", "20", "21", "22"], "params": {"maxdis": 0.1}},
}


def parse(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.strip().startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])
    return rows


def run_rolloff(binary: Path, cfg: dict, params: dict, work: Path, tag: str) -> Path:
    (work / "refs").write_text("\n".join(cfg["refs"]) + "\n")
    out = work / f"{tag}.out"
    par = {"genotypename": f"{cfg['prefix']}.geno", "snpname": f"{cfg['prefix']}.snp",
           "indivname": f"{cfg['prefix']}.ind", "poplistname": str(work / "refs"), "admixpop": cfg["admixpop"],
           "output": str(out), **params}
    if cfg["chromosomes"]:
        keep = set(cfg["chromosomes"])
        bad = [line.split()[0] for line in Path(cfg["prefix"] + ".snp").read_text().splitlines()
               if line.split()[1] not in keep]
        (work / "bad.snps").write_text("\n".join(bad) + "\n")
        par["badsnpname"] = str(work / "bad.snps")
    (work / f"{tag}.par").write_text("".join(f"{k}: {v}\n" for k, v in par.items()))
    log = subprocess.run([str(binary), "-p", str(work / f"{tag}.par")], capture_output=True, text=True,
                         check=True).stdout
    (work / f"{tag}.log").write_text(log)
    return out


def fit(path: Path) -> list[float | None]:
    out = subprocess.run(["/usr/bin/Rscript", str(HERE / "rexpfit_full.R"), str(path), "4", "0.5", "100", "1"],
                         capture_output=True, text=True, check=True).stdout.split()
    return [None if x == "NA" else float(x) for x in out]


def rexpfit_stock(path: Path, work: Path) -> float:
    """rexpfit.r itself (seeded the same way), its printed estimate of m."""
    script = work / "rexpfit_seeded.R"
    src = (BIN / "rexpfit.r").read_text().replace("library('DEoptim')", "library('DEoptim'); set.seed(1)")
    script.write_text('.libPaths(c("~/R/popmix-references", .libPaths()))\n' + src)
    subprocess.run(["/usr/bin/Rscript", str(script), str(path), str(work / "stockfit"), "4", "0.5", "100", "TRUE",
                    "FALSE", "FALSE"], capture_output=True, text=True, check=True)
    flog = (work / "stockfit.flog").read_text()
    return float(re.search(r"mean \(generations\):\s+(\S+)", flog).group(1))


def run(name: str) -> None:
    cfg = CONFIGS[name]
    params = {**PARAMS, **cfg.get("params", {})}
    print(f"== {name}", flush=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        hi = run_rolloff(BIN / "rolloff-hiprec", cfg, params, work, "hiprec")
        st = run_rolloff(BIN / "rolloff-stock", cfg, params, work, "stock")
        bins, stock_bins = parse(hi), parse(st)
        assert len(bins) == len(stock_bins)
        for h, s in zip(bins, stock_bins):  # noqa: B905
            assert all(abs(a - b) <= 5e-7 + 1e-12 * abs(a) for a, b in zip(h[1:4], s[1:4])), (h, s)  # noqa: B905
            assert h[4] == s[4]
        chroms = sorted({p.name.split(":")[-1] for p in work.glob("hiprec.out:*")}, key=int)
        jack = {c: parse(work / f"hiprec.out:{c}") for c in chroms}
        # chromosomes with no pairs have no data
        jack = {c: rows for c, rows in jack.items() if any(r[1] > 0 for r in rows)}
        full_fit = fit(hi)
        assert abs(rexpfit_stock(hi, work) - full_fit[2]) <= 5e-8 * abs(full_fit[2]) + 1e-12, "rexpfit.r disagrees"
        jack_fits = {c: fit(work / f"hiprec.out:{c}") for c in jack}
        snp_counts = Counter(line.split()[1] for line in Path(cfg["prefix"] + ".snp").read_text().splitlines())
        jin = work / "run.jin"
        jin.write_text("".join(f"{c} {snp_counts[c]} {jack_fits[c][2]!r}\n" for c in jack))
        subprocess.run([str(BIN / "dowtjack-hiprec"), "-i", str(jin), "-m", repr(full_fit[2]), "-o",
                        str(work / "run.jout")], capture_output=True, text=True, check=True)
        jout = (work / "run.jout").read_text()
        m = re.search(r"mean:\s+(\S+)\s+std error:\s+(\S+)", jout)
        n_used = re.search(r"snps processes \(outer \+ inner loop\) (\d+) (\d+)", (work / "hiprec.log").read_text())
    out = {
        "admixtools_version": "7.0.2",
        "prefix": os.path.relpath(cfg["prefix"], REPO) if cfg["prefix"].startswith(str(REPO)) else cfg["prefix"],
        "admixpop": cfg["admixpop"], "refs": cfg["refs"], "chromosomes": cfg["chromosomes"], "params": params,
        "bins": bins, "jackknife_bins": jack, "fit": full_fit, "jackknife_fits": jack_fits,
        "snp_weights": {c: snp_counts[c] for c in jack},
        "jackknife": [float(m.group(1)), float(m.group(2))],
        "pairs_visited": [int(n_used.group(1)), int(n_used.group(2))],
    }
    dest = HERE / "v7.0.2" / name
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "rolloff.json").write_text(json.dumps(out) + "\n")
    print(f"  {len(bins)} bins, decay {full_fit[2]:.4f}, jackknife {out['jackknife']}", flush=True)


if __name__ == "__main__":
    for name in sys.argv[1:] or list(CONFIGS):
        run(name)
