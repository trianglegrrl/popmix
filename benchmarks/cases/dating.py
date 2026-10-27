"""Benchmark cases: ROLLOFF and DATES."""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from ._common import CACHE, CLASSIC, HOME, Case, subset

# two chromosomes (the jackknife needs two), to 5 cM: rolloff takes about 7
# minutes for chromosome 22 alone, and over half an hour for 20-22 to 10 cM
ROLLOFF = {"admixed": "French", "refs": ["Sardinian", "Russia_Samara_EBA_Yamnaya"], "binsize": 0.001,
           "maxdis": 0.05, "chromosomes": ["21", "22"]}


def _rolloff_badsnps_path() -> Path:
    return CACHE / f"rolloff_chr{'_'.join(ROLLOFF['chromosomes'])}.badsnps"


def rolloff_badsnps() -> Path:
    """The subset's SNPs off ROLLOFF's chromosomes (rolloff's badsnpname, as
    generate_rolloff.py). Returns the subset."""
    src = subset()
    bad = _rolloff_badsnps_path()
    if not bad.exists():
        keep = set(ROLLOFF["chromosomes"])
        bad.write_text("".join(f"{s[0]}\n" for s in (line.split() for line in
                                                      Path(f"{src}.snp").read_text().splitlines())
                               if s[1] not in keep))
    return src


def _rolloff(i: Path, o: Path) -> list[str]:
    # tests/golden/rolloff/generate_rolloff.py's parameters
    (o / "refs").write_text("".join(f"{p}\n" for p in ROLLOFF["refs"]))
    par = o / "rolloff.par"
    par.write_text("".join(f"{k}: {v}\n" for k, v in {
        "genotypename": f"{i}.geno", "snpname": f"{i}.snp", "indivname": f"{i}.ind", "poplistname": o / "refs",
        "admixpop": ROLLOFF["admixed"], "output": o / "rolloff.out", "binsize": ROLLOFF["binsize"],
        "maxdis": ROLLOFF["maxdis"], "runmode": 1, "chithresh": 6.0, "zdipcorrmode": "YES", "checkmap": "YES",
        "seed": 1, "numchrom": 22, "jackknife": "YES", "badsnpname": _rolloff_badsnps_path()}.items()))
    return [str(CLASSIC / "rolloff-hiprec"), "-p", str(par)]


def _agree_rolloff(ours: Path, theirs: Path) -> dict:
    want = np.array([[float(x) for x in line.split()] for line in (theirs / "rolloff.out").read_text().splitlines()
                     if line.strip() and not line.lstrip().startswith("#")])
    import pandas as pd

    got = pd.read_csv(ours / "out.bins.tsv", sep="\t")
    # tests/golden/test_rolloff.py: curves rtol 1e-9, pair counts exact
    cols = {"zscore": 1, "regression": 2, "corr": 3}
    rel = {k: float(np.max(np.abs(got[k].to_numpy() - want[:, j]) / np.maximum(np.abs(want[:, j]), 1e-12)))
           for k, j in cols.items()} if len(got) == len(want) else {k: float("inf") for k in cols}
    counts = (len(got) == len(want) and bool((got["count"].to_numpy() == want[:, 4]).all())
              and np.allclose(got["d_cm"].to_numpy(), want[:, 0], rtol=1e-12, atol=0))
    return {"bins": len(want), **{f"{k}_rel": v for k, v in rel.items()}, "bins_equal": counts,
            "pass": counts and max(rel.values()) <= 1e-9}


ROLLOFF_CASE = Case(
    name="rolloff", analysis=f"ROLLOFF, {ROLLOFF['admixed']} = {' + '.join(ROLLOFF['refs'])}, "
                             f"to {ROLLOFF['maxdis'] * 100:g} cM, chromosomes {ROLLOFF['chromosomes'][0]}-"
                             f"{ROLLOFF['chromosomes'][-1]}, jackknife over them",
    reference="ADMIXTOOLS 7.0.2 rolloff (full-precision build)",
    popmix=lambda i, o, t: ["popmix", "--threads", str(t), "rolloff", "-i", str(i), "--admixed", ROLLOFF["admixed"],
                            "--refs", ",".join(ROLLOFF["refs"]), "--binsize", str(ROLLOFF["binsize"]),
                            "--maxdis", str(ROLLOFF["maxdis"]), "--chromosomes", ",".join(ROLLOFF["chromosomes"]),
                            "-o", str(o / "out")],
    ref=_rolloff, agree=_agree_rolloff, prepare=rolloff_badsnps,
)


DATES = Path(os.environ.get("POPMIX_DATES", HOME / "DATES"))


# Russia_Samara_EBA_Yamnaya mixes pseudo-haploid samples and a diploid genome
# (I6729.DG): --compat dates-v753 counts every sample diploid, as DATES does
# (DEVIATIONS.md, DATES-0004)
DATES_POPS = {"test": "French", "source1": "Sardinian", "source2": "Russia_Samara_EBA_Yamnaya"}


DATES_BINSIZE, DATES_MAXDIS, DATES_QBIN, DATES_LOVAL = 0.001, 0.5, 10, 0.45


# tests/golden/test_dates_v753.py: printed precision, plus dates_expfit's
# optimizer tolerance (gslprecision 1e-6 in r) for fitted generations
COV_TOL, CORR_TOL, GEN_TOL = 5e-10 + 1e-15, 5e-7 + 1e-12, 5e-4 + 1e-9


FIT_TOL = GEN_TOL + 1e-6 / DATES_BINSIZE


REPO = Path(__file__).resolve().parents[2]


DATES_BUILD = CACHE / "dates-v753-patched"


def dates_tools() -> Path:
    """DATES v753 built with tests/golden/dates/dates_numbins_fix.patch, as the
    goldens use, and its helper scripts on one PATH (DATES_BUILD/tools).

    The stock bin/dates reads an uninitialized buffer size and aborts at
    random ("Unable to allocate -1817970720 unit(s)"), in mesh mode too, and
    its dowtjack needs a libgsl this system no longer has. The build is from
    the DATES checkout's committed source, in a clean environment (conda's
    compiler flags break the link). Returns the subset."""
    tools = DATES_BUILD / "tools"
    if not (tools / "dowtjack").exists():
        src = DATES_BUILD / "src"
        DATES_BUILD.mkdir(parents=True, exist_ok=True)
        archive = subprocess.run(["git", "-C", str(DATES), "archive", "HEAD", "src"], check=True,
                                 capture_output=True).stdout
        subprocess.run(["tar", "-x", "-C", str(DATES_BUILD)], input=archive, check=True)
        # the archive carries prebuilt objects and libraries: build from source
        subprocess.run(["sh", "-c", "rm -f *.o nicksrc/*.o nicksrc/libnick.a nicklib dates dowtjack"], cwd=src,
                       check=True)
        subprocess.run(["patch", "-p1", "-i", str(REPO / "tests/golden/dates/dates_numbins_fix.patch")], cwd=src,
                       check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["make", "dates", "dowtjack"], cwd=src, check=True, capture_output=True,
                       env={"PATH": "/usr/local/bin:/usr/bin:/bin", "HOME": str(HOME)})
        tools.mkdir(exist_ok=True)
        for prog in [*(DATES / "bin").iterdir(), src / "dates", src / "dowtjack"]:
            link = tools / prog.name
            if link.is_symlink():
                link.unlink()
            link.symlink_to(prog)
    return subset()


def _dates(i: Path, o: Path) -> list[str]:
    p = DATES_POPS
    (o / "admixlist").write_text(f"{p['source1']} {p['source2']} {p['test']} {p['test']}\n")
    par = o / "dates.par"
    # tests/golden/dates/generate_dates_reference.sh's parameters
    par.write_text("".join(f"{k}: {v}\n" for k, v in {
        "genotypename": f"{i}.geno", "snpname": f"{i}.snp", "indivname": f"{i}.ind", "admixlist": "admixlist",
        "binsize": DATES_BINSIZE, "maxdis": DATES_MAXDIS, "seed": 77, "runmode": 1, "chithresh": 0.0,
        "mincount": 4, "qbin": DATES_QBIN, "runfit": "YES", "afffit": "YES", "lovalfit": DATES_LOVAL}.items()))
    tools = DATES_BUILD / "tools"
    # dates runs its helpers (dates_jackknife, dowtjack, ...) by name
    return ["env", f"PATH={tools}:{os.environ['PATH']}", str(tools / "dates"), "-p", par.name]


def _agree_dates(ours: Path, theirs: Path) -> dict:
    test = DATES_POPS["test"]
    res = theirs / test
    got = json.loads((ours / "dates.json").read_text())
    rows = [line.split() for line in (res / f"{test}.out").read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    # the JSON holds the first 20 bins
    ld = [(abs(b["n_pairs"] - int(r[4])), abs(b["weighted_ld"] - float(r[1])), abs(b["regression_coef"] - float(r[2])),
           abs(b["correlation"] - float(r[3]))) for b, r in zip(got["ld_bins"], rows)]  # noqa: B905
    full = float(re.search(r"mean \(generations\):\s+(\S+)", (res / f"expfit_{test}.flog").read_text()).group(1))
    blocks = sum(float(line.split()[1]) > 0 for line in (res / f"{test}.jin").read_text().splitlines())
    jk_est, jk_se = (float(x) for x in (res / f"{test}.jout").read_text().split())
    jk_tol = (blocks + 2) * FIT_TOL + GEN_TOL
    out = {"bins_compared": len(ld), "pairs_equal": all(d[0] == 0 for d in ld),
           "max_cov_abs": max(d[1] for d in ld), "max_corr_abs": max(max(d[2], d[3]) for d in ld),
           "fit_abs": abs(got["fit"]["decay_rate"] - full), "generations_abs": abs(got["generations"] - jk_est),
           "se_abs": abs(got["generations_se"] - jk_se)}
    out["pass"] = (out["pairs_equal"] and out["max_cov_abs"] <= COV_TOL and out["max_corr_abs"] <= CORR_TOL
                   and out["fit_abs"] <= FIT_TOL and out["generations_abs"] <= jk_tol and out["se_abs"] <= jk_tol)
    return out


def _popmix_dates(i: Path, o: Path, t: int) -> list[str]:
    p = DATES_POPS
    return ["sh", "-c", f"popmix --threads {t} dates -i {i} --test-pop {p['test']} --source1 {p['source1']} "
                        f"--source2 {p['source2']} --bin-width {DATES_BINSIZE * 100:g} "
                        f"--max-distance {DATES_MAXDIS * 100:g} --min-distance {DATES_LOVAL} --qbin {DATES_QBIN} "
                        f"--compat dates-v753 --format json > {o / 'dates.json'}"]


DATES_CASE = Case(
    name="dates", analysis=f"DATES, {DATES_POPS['test']} = {DATES_POPS['source1']} + {DATES_POPS['source2']}, "
                           f"FFT mesh (qbin {DATES_QBIN}), to {DATES_MAXDIS * 100:g} cM, jackknife",
    reference="DATES v753 (numbins fix, as the goldens)",
    popmix=_popmix_dates, ref=_dates, agree=_agree_dates, prepare=dates_tools,
)
