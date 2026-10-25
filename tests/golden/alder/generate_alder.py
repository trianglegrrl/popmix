#!/usr/bin/env python3
"""Generate the ALDER 1.03 reference for popmix's ALDER.

Usage: python tests/golden/alder/generate_alder.py [config...]

Runs ALDER 1.03 built to print full precision (alder-hiprec, from
build_alder.sh) and as released (alder-stock), checks that the two agree to
the stock build's printed precision, and writes
tests/golden/alder/v1.03/<config>/alder.json with, for each run of the
config (the main run, and a 1-reference run with each reference so that
every curve the admixture test uses is recorded):

  - params: the ALDER parameters beyond the file names
  - ld_corr_stops: the distance (cM) to which LD in the test population is
    correlated with LD in each reference, as ALDER summarizes it
  - fit_start: the distance (cM) from which the curve is fit
  - curve: the weighted LD curve (d in cM, weighted LD, bin count), the last
    row the inter-chromosome affine term (d = inf)
  - jackknife_curves: the same with each chromosome left out
  - fits: for each curve the run computes (keyed "2-ref A B", "1-ref A"),
    for each fit start ALDER tries, decay, amp_tot, amp_exp and amp_aff,
    each as the full-data value, the jackknife standard error and the
    jackknife replicates
  - test: the admixture test (2-reference runs): z-scores of both 1-reference
    curves and of the 2-reference curve, decay differences, p-value, and
    whether it succeeds
  - mixture_fraction: the 1-reference mixture-fraction lower bound
"""

from __future__ import annotations

import glob
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AADR = "/references/aadr/v62.0_1240k_public"
SMALL = str(REPO / "tests/golden/dates/family_chr20_22")
BIN = Path(os.environ.get("ALDER_DIR", Path.home() / ".local/share/popmix-references/alder-1.03"))

CONFIGS = {
    "family": {"prefix": SMALL, "admixpop": "Simulation.1", "refpops": ["French", "Yoruba"]},
    "aadr_asw": {"prefix": AADR, "admixpop": "ASW.DG", "refpops": ["YRI.DG", "CEU.DG"]},
    # a pseudo-haploid ancient reference with no data at some SNPs (ALDER-0001)
    "aadr_asw_lbk": {"prefix": AADR, "admixpop": "ASW.DG", "refpops": ["YRI.DG", "Germany_EN_LBK.AG"]},
    "aadr_asw_mindis": {"prefix": AADR, "admixpop": "ASW.DG", "refpops": ["YRI.DG", "CEU.DG"],
                            "params": {"mindis": 0.008}},
}


def num(x: str) -> float:
    return math.inf if x in ("inf", "INFINITY") else (math.nan if x == "nan" else float(x))


def parse_fits(log: str) -> list[dict]:
    fits = []
    for block in re.split(r"(?=---- fit on data from )", log)[1:]:
        head = re.match(r"---- fit on data from (\S+) to (\S+) cM \((.*?)\) ----", block)
        fit = {"mindis_cM": float(head.group(1)), "maxdis_cM": float(head.group(2)),
               "inter_chrom_affine": head.group(3).startswith("using")}
        for var in ("decay", "amp_tot", "amp_exp", "amp_aff"):
            m = re.search(rf"d>\S+\s+{var}:(\S+)(?: \+/- (\S+))?", block)
            if not m:
                break
            reps = re.search(rf"d>\S+\s+{var} jackknife reps\t(.*)", block)
            fit[var] = {"value": num(m.group(1)), "se": None if m.group(2) is None else num(m.group(2)),
                        "reps": [num(x) for x in reps.group(1).split()] if reps else None}
        if "decay" in fit:
            fits.append(fit)
    # the comparison section repeats the headers without fits
    return fits


def parse_curve(path: str) -> list[list[float]]:
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("# last row") or line.startswith("#   d"):
            continue
        f = line.lstrip("# ").split()
        rows.append([num(f[0]), num(f[1]), int(f[2])])
    return rows


def run(binary: Path, cfg: dict, refpops: list[str], work: Path, tag: str) -> dict:
    par = {
        "genotypename": f"{cfg['prefix']}.geno", "snpname": f"{cfg['prefix']}.snp",
        "indivname": f"{cfg['prefix']}.ind", "admixpop": cfg["admixpop"], "refpops": ";".join(refpops),
        "raw_outname": str(work / f"{tag}.raw"), "jackknife": "YES", "print_jackknife_fits": "YES",
        **{k: str(v) for k, v in cfg.get("params", {}).items()},
    }
    (work / f"{tag}.par").write_text("".join(f"{k}: {v}\n" for k, v in par.items()))
    log = subprocess.run([str(binary), "-p", str(work / f"{tag}.par")], capture_output=True, text=True,
                         check=True).stdout
    (work / f"{tag}.log").write_text(log)
    stops = dict(re.findall(r"^\s*(\S+) (\S+)(?: --> .*)?$",
                            log.split("Summary of LD correlation results")[1].split("==>")[0], re.M)) \
        if "Summary of LD correlation results" in log else {}
    stops = {k: num(v) for k, v in stops.items() if k in refpops}
    fit_start = re.search(r"Weighted LD curves are fit starting at (\S+) cM", log)
    out = {
        "refpops": refpops,
        "ld_corr_stops": stops,
        "fit_start": None if fit_start is None else float(fit_start.group(1)),
        "curve": parse_curve(str(work / f"{tag}.raw")),
        "jackknife_curves": {Path(p).name.split("-")[-1]: parse_curve(p)
                             for p in sorted(glob.glob(str(work / f"{tag}.raw-*")))},
        # fits of each curve computed in the run, keyed by its weights
        # ("2-ref A B", "1-ref A")
        "fits": {f"{n}-ref {w.strip()}": parse_fits(body.split("Comparing curves to test")[0])
                 for n, w, body in re.findall(r"\*\*\* Computing (\d)-ref weighted LD with weights ([^*]+) \*\*\*"
                                             r"(.*?)(?=\*\*\* Computing |\Z)", log, re.S)},
    }
    if "Admixture test summary" in log:
        t = log.split("Admixture test summary")[1]
        z = lambda label: num(re.search(label + r"\s+(\S+)", t).group(1))  # noqa: E731
        result = re.search(r"Test (SUCCEEDS|FAILS) \(z=(\S+), p=(\S+)\)", t)
        oneref = re.findall(r"1-ref decay z-score:\s+(\S+)\n\s+1-ref amp_exp z-score:\s+(\S+)", t)
        diffs = re.findall(r"^\s+(1-ref .*?) z-score:\s+(\S+)\s+\((\S+)%\)", t, re.M)
        out["test"] = {
            "oneref_z": [[num(a), num(b)] for a, b in oneref],
            "tworef_z": [z("2-ref decay z-score:"), z("2-ref amp_exp z-score:")],
            "decay_diffs": [{"label": lab, "z": num(a), "percent": num(b)} for lab, a, b in diffs],
            "consistent": "YES: decay rates are consistent" in t,
            "succeeds": result.group(1) == "SUCCEEDS",
            "z": num(result.group(2)),
            "p": num(result.group(3)),
        }
    mix = re.search(r"Mixture fraction % lower bound \(assuming admixture\): (\S+) \+/- (\S+)", log)
    if mix:
        out["mixture_fraction"] = [num(mix.group(1)), num(mix.group(2))]
    return out


def close(a, b, digits: int) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        if math.isinf(a) or math.isinf(b) or math.isnan(a) or math.isnan(b):
            return (math.isinf(a) and math.isinf(b)) or (math.isnan(a) and math.isnan(b))
        return abs(a - b) <= 0.5 * 10 ** -digits + 1e-12 * abs(a)
    return True


def check_stock(hiprec: dict, stock: dict) -> None:
    """The stock build prints the curve with 8 decimals and decay with 2."""
    assert len(hiprec["curve"]) == len(stock["curve"])
    for h, s in zip(hiprec["curve"], stock["curve"]):  # noqa: B905
        assert h[2] == s[2] and close(h[1], s[1], 8), (h, s)
    assert hiprec["fits"].keys() == stock["fits"].keys()
    for key in hiprec["fits"]:
        for fh, fs in zip(hiprec["fits"][key], stock["fits"][key]):  # noqa: B905
            assert close(fh["decay"]["value"], fs["decay"]["value"], 2), (fh["decay"], fs["decay"])
            assert close(fh["amp_exp"]["value"], fs["amp_exp"]["value"], 8)
    assert hiprec.get("test", {}).get("succeeds") == stock.get("test", {}).get("succeeds")


def main(names: list[str]) -> None:
    for name in names or list(CONFIGS):
        cfg = CONFIGS[name]
        print(f"== {name}", flush=True)
        runs = {}
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for key, refs in [("main", cfg["refpops"])] + [(f"oneref_{r}", [r]) for r in cfg["refpops"]]:
                hiprec = run(BIN / "alder-hiprec", cfg, refs, work, key + "_h")
                stock = run(BIN / "alder-stock", cfg, refs, work, key + "_s")
                check_stock(hiprec, stock)
                runs[key] = hiprec
        out = {"alder_version": "1.03", "prefix": os.path.relpath(cfg["prefix"], REPO)
               if cfg["prefix"].startswith(str(REPO)) else cfg["prefix"], "admixpop": cfg["admixpop"],
               "refpops": cfg["refpops"], "params": cfg.get("params", {}), "runs": runs}
        dest = HERE / "v1.03" / name
        dest.mkdir(parents=True, exist_ok=True)
        (dest / "alder.json").write_text(json.dumps(out) + "\n")
        main_run = runs["main"]
        print(f"  fit start {main_run['fit_start']} cM, {len(main_run['curve'])} bins, "
              f"test {main_run.get('test', {}).get('succeeds')}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
