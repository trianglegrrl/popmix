#!/usr/bin/env python3
"""Scaling curves (TASKS.md E5): time and peak memory against the size of
the input, popmix and the reference tool on the same inputs.

Curves, each on inputs cut from AADR v66.p1 by AdmixTools 8 convertf:

- populations: f2 extraction for 4 to 497 populations of 10 samples each
  (every group of AADR with 10 or more samples);
- snps: f2 extraction for the 16 benchmark populations on the first
  1/16 to all of the 1,233,013 SNPs;
- samples: PCA and KING kinship on 250 to 3,970 present-day samples (all of
  AADR's);
- threads: popmix at 1, 2, 4, 8 and 16 threads (and PLINK 2 for kinship);
- full_aadr: f2 extraction of the 16 populations from all of AADR (23,089
  samples), and popmix from AADR's own TGENO file.

Every point checks agreement as `run.py` does. Points already in the output
file are skipped, so an interrupted run resumes.

Usage: python benchmarks/scaling.py [--curves populations,snps] [--repeats 3]
                                    [--threads 16] [--out scaling.json]
"""

from __future__ import annotations

import argparse
import datetime
import json
import random
import shutil
import socket
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

sys.path.insert(0, str(Path(__file__).parent))

from cases import (
    _R_F2,  # noqa: E402
    AADR,
    CACHE,
    CONVERTF,
    PLINK2,
    POPS,
    SMARTPCA,
    _agree_f2,
    _agree_king,
    _agree_pca,
    _r,
    subset,
)
from harness import machine, measure  # noqa: E402
from run import _plain, versions  # noqa: E402

SCALING = CACHE / "scaling"
ANNO = AADR.parent / "v66.p1_1240K.aadr.PUB.anno"
PER_POP = 10
POP_COUNTS = (4, 8, 16, 32, 64, 128, 256, 497)
SNP_FRACTIONS = (16, 8, 4, 2, 1)  # 1/k of the SNPs
SAMPLE_COUNTS = (250, 500, 1000, 2000, 3970)
THREADS = (1, 2, 4, 8, 16)
THREAD_POPS, THREAD_SAMPLES = 128, 2000
# admixtools2 above this many populations takes hours; popmix alone beyond it
REF_MAX_POPS = 256
SEED = 1


def _ind_rows() -> list[tuple[str, str, str]]:
    return [tuple(line.split()) for line in Path(f"{AADR}.ind").read_text().splitlines()]


def pop_groups() -> list[str]:
    """Groups of AADR with PER_POP or more samples, in a fixed shuffled order,
    so the first n of them are a nested series."""
    counts: dict[str, int] = {}
    for _, _, group in _ind_rows():
        counts[group] = counts.get(group, 0) + 1
    groups = sorted(g for g, c in counts.items() if c >= PER_POP and "Ignore" not in g)
    random.Random(SEED).shuffle(groups)
    return groups


def present_day_samples() -> list[str]:
    """AADR's present-day samples (date 0 BP), in a fixed shuffled order."""
    rows = [line.split("\t") for line in ANNO.read_text().splitlines()[1:]]
    ids = sorted(r[0] for r in rows if r[10].strip() == "0")
    random.Random(SEED).shuffle(ids)
    return ids


def _convertf(name: str, keep: Callable[[str, str], bool], fmt: str = "PACKEDANCESTRYMAP",
              src: Path = AADR, label: Callable[[str], str] = lambda g: g, **extra: str) -> Path:
    """convertf's copy of src with only the samples `keep(id, group)` accepts,
    labelled `label(group)` (the others relabelled Ignore). Made once, in the
    cache."""
    out = SCALING / name
    done = Path(f"{out}.bed" if fmt == "PACKEDPED" else f"{out}.geno")
    if done.exists():
        return out
    SCALING.mkdir(parents=True, exist_ok=True)
    rows = [line.split() for line in Path(f"{src}.ind").read_text().splitlines()]
    ind = SCALING / f"{name}.in.ind"
    ind.write_text("".join(f"{i} {s} {label(g) if keep(i, g) else 'Ignore'}\n" for i, s, g in rows))
    ext = ("bed", "bim", "fam") if fmt == "PACKEDPED" else ("geno", "snp", "ind")
    keys = {"genotypename": f"{src}.geno", "snpname": f"{src}.snp", "indivname": ind, "outputformat": fmt,
            "genotypeoutname": f"{out}.{ext[0]}", "snpoutname": f"{out}.{ext[1]}",
            "indivoutname": f"{out}.{ext[2]}", **extra}
    par = SCALING / f"{name}.par"
    par.write_text("".join(f"{k}: {v}\n" for k, v in keys.items()))
    subprocess.run([CONVERTF, "-p", str(par)], cwd=SCALING, check=True, stdout=subprocess.DEVNULL)
    return out


@dataclass
class Point:
    curve: str
    x: int
    analysis: str
    prepare: Callable[[], Path]
    popmix: Callable[[Path, Path, int], list[str]]
    ref: Callable[[Path, Path, int], list[str]] | None
    agree: Callable[[Path, Path], dict] | None
    threads: tuple[int, ...] | None = None  # None: --threads
    ref_threads: tuple[int, ...] = (0,)  # 0: the reference's own default

    @property
    def key(self) -> str:
        return f"{self.analysis}:{self.x}"


def _f2_point(curve: str, x: int, pops: list[str], prepare: Callable[[], Path], with_ref: bool = True,
              threads: tuple[int, ...] | None = None) -> Point:
    return Point(
        curve, x, "extract_f2", prepare,
        popmix=lambda i, o, t: ["popmix", "--threads", str(t), "extract-f2", "-i", str(i), "-o", str(o / "f2"),
                                "--populations", ",".join(pops), "--pseudohaploid", "admixtools2"],
        ref=(lambda i, o, _t: _r(_R_F2, o, str(i), ",".join(pops))) if with_ref else None,
        agree=_agree_f2 if with_ref else None, threads=threads)


def _pca_point(curve: str, x: int, prepare: Callable[[], Path], threads: tuple[int, ...] | None = None) -> Point:
    def smartpca(i: Path, o: Path, _t: int) -> list[str]:
        par = o / "smartpca.par"
        par.write_text("".join(f"{k}: {v}\n" for k, v in {
            "genotypename": f"{i}.geno", "snpname": f"{i}.snp", "indivname": f"{i}.ind",
            "evecoutname": o / "out.evec", "evaloutname": o / "out.eval", "numoutevec": 10,
            "numoutlieriter": 0}.items()))
        return [SMARTPCA, "-p", str(par)]

    return Point(curve, x, "pca", prepare,
                 popmix=lambda i, o, t: ["popmix", "--threads", str(t), "pca", "-i", str(i), "--components", "10",
                                         "--outlier-iterations", "0", "-o", str(o / "out")],
                 ref=smartpca, agree=_agree_pca, threads=threads)


def _king_point(curve: str, x: int, prepare: Callable[[], Path], threads: tuple[int, ...] | None = None) -> Point:
    return Point(curve, x, "kinship", prepare,
                 popmix=lambda i, o, t: ["popmix", "--threads", str(t), "kinship", "-i", str(i), "--method", "king",
                                         "-o", str(o / "kin.tsv")],
                 ref=lambda i, o, t: [PLINK2, "--bfile", str(i), "--make-king-table", "--make-king", "triangle",
                                      "bin", "--threads", str(t or 16), "--out", str(o / "king")],
                 agree=_agree_king, threads=threads, ref_threads=threads or (0,))


def _pops_input(n: int) -> Callable[[], Path]:
    def prepare() -> Path:
        chosen = set(pop_groups()[:n])
        seen: dict[str, int] = {}

        def keep(_i: str, g: str) -> bool:
            if g not in chosen or seen.get(g, 0) >= PER_POP:
                return False
            seen[g] = seen.get(g, 0) + 1
            return True

        return _convertf(f"pops{n}", keep)
    return prepare


# The sample series as one population: smartpca stops at more than 100
# populations, and neither PCA nor KING uses them
PRESENT_DAY_LABEL = "PresentDay"


def _samples_input(n: int, fmt: str = "PACKEDANCESTRYMAP") -> Callable[[], Path]:
    def prepare() -> Path:
        chosen = set(present_day_samples()[:n])
        name = f"samples{n}_1pop" + ("_bed" if fmt == "PACKEDPED" else "")
        return _convertf(name, lambda i, _g: i in chosen, fmt, label=lambda _g: PRESENT_DAY_LABEL)
    return prepare


def _snps_input(k: int) -> Callable[[], Path]:
    def prepare() -> Path:
        src = subset()
        if k == 1:
            return src
        snps = [line.split()[0] for line in Path(f"{src}.snp").read_text().splitlines()]
        bad = SCALING / f"snps{k}.bad"
        SCALING.mkdir(parents=True, exist_ok=True)
        bad.write_text("".join(f"{s}\n" for s in snps[len(snps) // k:]))
        return _convertf(f"snps{k}", lambda _i, _g: True, src=src, badsnpname=str(bad))
    return prepare


def _full_packed() -> Path:
    return _convertf("full", lambda _i, _g: True)


def points() -> list[Point]:
    groups = pop_groups()
    out = [_f2_point("populations", n, groups[:n], _pops_input(n), with_ref=n <= REF_MAX_POPS)
           for n in POP_COUNTS]
    out += [_f2_point("snps", k, POPS, _snps_input(k)) for k in SNP_FRACTIONS]
    out += [_pca_point("samples", n, _samples_input(n)) for n in SAMPLE_COUNTS]
    out += [_king_point("samples", n, _samples_input(n, "PACKEDPED")) for n in SAMPLE_COUNTS]
    out += [_f2_point("threads", THREAD_POPS, groups[:THREAD_POPS], _pops_input(THREAD_POPS), threads=THREADS),
            _pca_point("threads", THREAD_SAMPLES, _samples_input(THREAD_SAMPLES), threads=THREADS),
            _king_point("threads", THREAD_SAMPLES, _samples_input(THREAD_SAMPLES, "PACKEDPED"), threads=THREADS)]
    out += [_f2_point("full_aadr", 23089, POPS, _full_packed),
            Point("full_aadr", 23089, "extract_f2_tgeno", lambda: AADR,
                  popmix=lambda i, o, t: ["popmix", "--threads", str(t), "extract-f2", "-i", str(i), "-o",
                                          str(o / "f2"), "--populations", ",".join(POPS)],
                  ref=None, agree=None)]
    return out


def _size(prefix: Path) -> dict:
    """Samples, SNPs and populations of a prefix (EIGENSTRAT or PLINK)."""
    ind = Path(f"{prefix}.ind") if Path(f"{prefix}.ind").exists() else Path(f"{prefix}.fam")
    snp = Path(f"{prefix}.snp") if Path(f"{prefix}.snp").exists() else Path(f"{prefix}.bim")
    rows = [line.split() for line in ind.read_text().splitlines()]
    col = 2 if ind.suffix == ".ind" else 0
    used = [r for r in rows if r[col] != "Ignore"]
    with open(snp) as f:
        n_snps = sum(1 for _ in f)
    return {"samples": len(used), "snps": n_snps, "populations": len({r[col] for r in used})}


def _fresh(directory: Path) -> Path:
    shutil.rmtree(directory, ignore_errors=True)
    directory.mkdir(parents=True)
    return directory


def reused_reference(earlier: dict, curve: str, key: str, work: Path) -> dict:
    """The reference timings of a point from an earlier report, whose outputs
    are still in `work`, each marked with the earlier date."""
    point = earlier.get("curves", {}).get(curve, {}).get(key)
    if point is None:
        raise SystemExit(f"the earlier report has no point {curve}/{key}")
    out = {}
    for t, runs in point["reference"].items():
        ref_dir = work / f"reference-t{t}"
        if not ref_dir.is_dir() or not any(ref_dir.iterdir()):
            raise SystemExit(f"{curve}/{key}: no reference outputs in {ref_dir}")
        out[t] = {**runs, "reused_from": earlier["date"]}
    return out


def run_point(p: Point, threads: int, repeats: int, earlier: dict | None = None) -> dict:
    source = p.prepare()
    work = SCALING / "runs" / p.curve / p.key.replace(":", "-")
    result: dict = {"x": p.x, "analysis": p.analysis, "input": _size(source), "reference": {}, "popmix": {}}
    ref_dir = None
    if p.ref is not None and earlier is not None:
        result["reference"] = reused_reference(earlier, p.curve, p.key, work)
        ref_dir = work / f"reference-t{list(result['reference'])[-1]}"
    elif p.ref is not None:
        for t in p.ref_threads:
            ref_dir = _fresh(work / f"reference-t{t}")
            print(f"  {p.key}: reference{f', {t} threads' if t else ''}", flush=True)
            result["reference"][str(t)] = measure(p.ref(source, ref_dir, t), ref_dir, repeats,
                                                  log=ref_dir / "run.log").summary()
    for t in p.threads or (threads,):
        ours = _fresh(work / f"popmix-t{t}")
        print(f"  {p.key}: popmix, {t} threads", flush=True)
        m = measure(p.popmix(source, ours, t), ours, repeats, log=ours / "run.log")
        result["popmix"][str(t)] = m.summary()
        if p.agree is not None and ref_dir is not None:
            result["popmix"][str(t)]["agreement"] = p.agree(ours, ref_dir)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--curves", default="populations,snps,samples,threads,full_aadr")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--threads", type=int, default=16, help="popmix threads on the size curves")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--reuse-reference", type=Path, metavar="REPORT",
                        help="take the reference timings from this earlier report and its outputs from the "
                             "run directories, and run only popmix")
    args = parser.parse_args()
    earlier = json.loads(args.reuse_reference.read_text()) if args.reuse_reference else None
    today = datetime.date.today().isoformat()
    out = args.out or Path(__file__).parent / "results" / f"scaling-{socket.gethostname()}-{today}.json"
    report = json.loads(out.read_text()) if out.exists() else {
        "date": today, "machine": machine(), "versions": versions(), "repeats": args.repeats, "curves": {}}
    curves = args.curves.split(",")
    for p in points():
        if p.curve not in curves:
            continue
        done = report["curves"].setdefault(p.curve, {})
        if p.key in done:
            continue
        print(f"{p.curve}: {p.key}", flush=True)
        done[p.key] = run_point(p, args.threads, args.repeats, earlier)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=1, default=_plain) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
