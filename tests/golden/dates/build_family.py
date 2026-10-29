#!/usr/bin/env python3
"""Build the small family data from AADR v66.p1 and a simulated admixture.

Writes, in PACKEDANCESTRYMAP (by AdmixTools 8 convertf):

    tests/golden/dates/family_chr20_22.{geno,snp,ind}            (in the repository)
    tests/golden/dates/example_data/family_packed.{geno,snp,ind}  (every autosome; gitignored)
    tests/golden/dates/family_truth.json                          (the simulation and its truth)

Both hold the same 50 samples:

- French (24) and Yoruba (21): HGDP genomes of AADR v66.p1 (CC0,
  doi:10.7910/DVN/FFIDCW) labelled French.DG and Yoruba.DG (diploid
  genotypes from whole-genome sequence), genotypes unchanged;
- Simulation (NA00000) and Simulation.1 (NA00001-NA00004): admixed samples
  simulated here, 20% French and 80% Yoruba, in one pulse 50 generations
  ago.

and the same SNPs: every third autosomal SNP of the 1240K panel (rows 0, 3,
6, ... of its .snp), which keeps the data the size of what they replace;
family_chr20_22 is family_packed's chromosomes 20-22.

The simulation. Each simulated sample is two haplotypes, each made
chromosome by chromosome. Ancestry changes along a haplotype at the
crossovers since the admixture: a Poisson process of rate 50 per Morgan on
the AADR genetic map (the .snp file's third column, in Morgans), from the
chromosome's first SNP to its last. Each segment between crossovers is
French with probability 0.2, else Yoruba, independently, so ancestry LD
decays as exp(-50 d), the model DATES and ROLLOFF fit. A segment's alleles
come from one source genome of its ancestry, chosen at random, never a
reference sample:

- French: the last 4 HGDP French genomes by sample ID (HGDP00536-539);
- Yoruba: 30 genomes of AADR's 1000 Genomes Yoruba (YRI.DG, Yoruba in
  Ibadan), in sample-ID order, each unrelated (KING-robust kinship below
  2^-4.5, third degree) to every HGDP Yoruba and to the sources before it.

A pool of sources this large keeps the simulated population from looking
drifted: with a handful of source genomes, their own sampling noise is
shared by every simulated sample and outweighs the admixture's negative
f3(Simulation.1; French, Yoruba). v66.p1's other French.DG and Yoruba.DG
genomes (SGDP and Meyer et al. 2012) are HGDP individuals sequenced again:
their calls agree with one HGDP genome's at over 99.9% of SNPs, so they are
left out.

The sources are unphased, so at each SNP the haplotype takes one of the
source's two alleles at random (a missing call stays missing); the sample's
genotype is the sum of its two haplotypes, missing if either is. numpy's
default_rng(SEED) draws everything, in a fixed order (sample, haplotype,
chromosome).

The truth (family_truth.json) records the parameters, the samples, and
per haplotype the number of tracts and the French share of its genetic
length, genome-wide and on chromosomes 20-22.

Input: the AADR v66.p1 subset tests/golden/aadr/build_aadr_v66p1.py writes
in the reference cache (~/.cache/popmix-references/aadr_v66p1/aadr_v66p1),
checked against its manifest.json, read as PACKEDANCESTRYMAP (the format
of the CONVERTF documentation: a 48-byte-or-longer header record, then one
record per SNP of 2 bits per sample, first sample in the high bits, 3 for
missing).

Usage: python tests/golden/dates/build_family.py [output directory]

(default: tests/golden/dates; tests/golden/test_family_build.py builds
into a temporary directory and checks the result is byte-identical)
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CACHE = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1"
AADR = CACHE / "aadr_v66p1"
CONVERTF = os.environ.get("POPMIX_CONVERTF",
                          str(Path.home() / ".local/share/popmix-references/admixtools-8.0.2/bin/convertf"))

SEED = 20261029
GENERATIONS = 50
FRACTION = {"French": 0.2, "Yoruba": 0.8}
THIN = 3
# AADR label -> label here
LABELS = {"French.DG": "French", "Yoruba.DG": "Yoruba"}
YRI = "YRI.DG"
SIMULATED = [("NA00000", "Simulation")] + [(f"NA0000{i}", "Simulation.1") for i in range(1, 5)]
# HGDP French genomes, the last by sample ID, whose alleles the simulated
# samples take; the others are the reference samples
N_FRENCH_SOURCES = 4
# YRI genomes whose alleles the simulated samples' Yoruba segments take
N_YRI_SOURCES = 30
# KING-robust kinship of third-degree relatives and closer (Manichaikul et al. 2010)
RELATED = 2 ** -4.5
# Genome pairs whose calls agree at more than this share of the SNPs both
# call are one person sequenced twice (unrelated pairs agree at about 60%).
DUPLICATE_CONCORDANCE = 0.95


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while chunk := f.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def read_aadr() -> tuple[list[tuple[str, str, str]], list[str], np.ndarray]:
    """The French.DG, Yoruba.DG and YRI.DG samples, every SNP line, genotypes SNPs x samples (9 missing)."""
    manifest = json.loads((CACHE / "manifest.json").read_text(encoding="utf-8"))
    for ext in ("geno", "snp", "ind"):
        if md5(Path(f"{AADR}.{ext}")) != manifest["outputs"][ext]["md5"]:
            raise SystemExit(f"{AADR}.{ext} differs from its manifest: rerun build_aadr_v66p1.py")
    inds = [tuple(line.split()) for line in Path(f"{AADR}.ind").read_text(encoding="utf-8").splitlines()]
    snps = Path(f"{AADR}.snp").read_text(encoding="utf-8").splitlines()
    keep = [k for k, (_, _, pop) in enumerate(inds) if pop in LABELS or pop == YRI]
    rlen = max(48, (len(inds) + 3) // 4)
    raw = np.fromfile(f"{AADR}.geno", dtype=np.uint8)
    header = raw[:rlen].tobytes().split(b"\0")[0].split()
    if header[0] != b"GENO" or int(header[1]) != len(inds) or int(header[2]) != len(snps):
        raise SystemExit(f"{AADR}.geno: unexpected header {header}")
    rec = raw[rlen:].reshape(len(snps), rlen)
    cols = np.array(keep)
    geno = (rec[:, cols // 4] >> (6 - 2 * (cols % 4)).astype(np.uint8)) & 3
    geno[geno == 3] = 9
    return [inds[k] for k in keep], snps, geno


def concordance(a: np.ndarray, b: np.ndarray) -> float:
    both = (a != 9) & (b != 9)
    return float(np.mean(a[both] == b[both]))


def king(geno: np.ndarray) -> np.ndarray:
    """KING-robust kinship of every pair of columns (Manichaikul et al. 2010), on SNPs both call."""
    called = (geno != 9).astype(np.float32)
    het = (geno == 1).astype(np.float32)
    hom0 = (geno == 0).astype(np.float32)
    hom2 = (geno == 2).astype(np.float32)
    hethet = het.T @ het
    ibs0 = hom0.T @ hom2 + hom2.T @ hom0
    het_i = het.T @ called  # [i, j]: heterozygous in i, called in j
    het_j = het_i.T
    low = np.minimum(het_i, het_j)
    return (hethet - 2 * ibs0) / (2 * low) + 0.5 - (het_i + het_j) / (4 * low)


def largest_kinship(geno: np.ndarray) -> float:
    kin = king(geno)
    return float(kin[~np.eye(len(kin), dtype=bool)].max())


def choose_samples(inds, geno) -> tuple[dict, dict, list[str], dict]:
    """Reference samples (HGDP), sources, the duplicates of HGDP genomes, and the largest kinships."""
    index = {s: k for k, (s, _, _) in enumerate(inds)}
    hgdp = {label: sorted(s for s, _, q in inds if q == pop and s.startswith("HGDP")) for pop, label in LABELS.items()}
    sources = {"French": hgdp["French"][-N_FRENCH_SOURCES:]}
    reference = {"French": [s for s, _, _ in inds if s in hgdp["French"][:-N_FRENCH_SOURCES]],
                 "Yoruba": [s for s, _, _ in inds if s in hgdp["Yoruba"]]}
    duplicates = []
    for pop, label in LABELS.items():
        for s, _, q in inds:
            if q != pop or s in hgdp[label]:
                continue
            twins = [r for r in hgdp[label]
                     if concordance(geno[:, index[s]], geno[:, index[r]]) > DUPLICATE_CONCORDANCE]
            if len(twins) != 1:
                raise SystemExit(f"{s} is not one HGDP genome sequenced again: {twins}")
            duplicates.append(f"{s} = {twins[0]}")
    yri = sorted(s for s, _, q in inds if q == YRI)
    names = reference["Yoruba"] + yri
    kin = king(geno[:, [index[s] for s in names]])
    pos = {s: k for k, s in enumerate(names)}
    kept: list[str] = []
    for s in yri:
        if max(kin[pos[s], pos[o]] for o in reference["Yoruba"] + kept) < RELATED:
            kept.append(s)
        if len(kept) == N_YRI_SOURCES:
            break
    sources["Yoruba"] = kept
    largest = {p: largest_kinship(geno[:, [index[s] for s in reference[p] + sources[p]]]) for p in reference}
    return reference, sources, duplicates, largest


def simulate(chrom: np.ndarray, morgans: np.ndarray, src: dict[str, np.ndarray], rng: np.random.Generator):
    """One haplotype: alleles (0/1, 9 missing) and per chromosome (tracts, French length, length)."""
    hap = np.full(len(chrom), 9, dtype=np.uint8)
    stats = {}
    for c in np.unique(chrom):
        rows = np.flatnonzero(chrom == c)
        pos = morgans[rows]
        lo, hi = pos[0], pos[-1]
        breaks = np.sort(rng.uniform(lo, hi, rng.poisson(GENERATIONS * (hi - lo))))
        edges = np.concatenate([[lo], breaks, [hi]])
        french = rng.random(len(edges) - 1) < FRACTION["French"]
        which = np.where(french, rng.integers(0, src["French"].shape[1], len(french)),
                         rng.integers(0, src["Yoruba"].shape[1], len(french)))
        seg = np.clip(np.searchsorted(breaks, pos, side="right"), 0, len(french) - 1)
        u = rng.random(len(rows))
        for k in range(len(french)):
            at = seg == k
            g = (src["French"] if french[k] else src["Yoruba"])[rows[at], which[k]]
            hap[rows[at]] = np.where(g == 9, 9, (u[at] < g / 2.0).astype(np.uint8))
        stats[int(c)] = (len(french), float(np.sum(np.diff(edges)[french])), float(hi - lo))
    return hap, stats


def write_packed(prefix: Path, inds: list[tuple[str, str, str]], snps: list[str], geno: np.ndarray, work: Path):
    """EIGENSTRAT text, converted to PACKEDANCESTRYMAP by convertf."""
    tmp = work / prefix.name
    Path(f"{tmp}.ind").write_text("".join(f"{s} {x} {p}\n" for s, x, p in inds), encoding="utf-8")
    Path(f"{tmp}.snp").write_text("".join(f"{line}\n" for line in snps), encoding="utf-8")
    text = (geno + ord("0")).astype(np.uint8)
    text[geno == 9] = ord("9")
    lines = np.concatenate([text, np.full((len(snps), 1), ord("\n"), dtype=np.uint8)], axis=1)
    lines.tofile(f"{tmp}.geno")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    par = work / f"{prefix.name}.par"
    par.write_text(f"genotypename: {tmp}.geno\nsnpname: {tmp}.snp\nindivname: {tmp}.ind\n"
                   f"outputformat: PACKEDANCESTRYMAP\ngenooutfilename: {prefix}.geno\n"
                   f"snpoutfilename: {prefix}.snp\nindoutfilename: {prefix}.ind\n", encoding="utf-8")
    # convertf writes a profile (gmon.out) where it runs
    subprocess.run([CONVERTF, "-p", str(par)], check=True, stdout=subprocess.DEVNULL, cwd=work)


def main(out: Path = HERE) -> None:
    small, full, truth_path = out / "family_chr20_22", out / "example_data" / "family_packed", out / "family_truth.json"
    inds, snps, geno = read_aadr()
    chrom_all = np.array([int(line.split()[1]) for line in snps])
    rows = np.flatnonzero((np.arange(len(snps)) % THIN == 0) & (chrom_all <= 22))
    snps = [snps[r] for r in rows]
    geno = geno[rows]
    chrom = chrom_all[rows]
    reference, sources, duplicates, largest = choose_samples(inds, geno)
    morgans = np.array([float(line.split()[2]) for line in snps])
    index = {s: k for k, (s, _, _) in enumerate(inds)}
    src = {p: geno[:, [index[s] for s in sources[p]]] for p in sources}

    rng = np.random.default_rng(SEED)
    sim_geno, samples = [], []
    for sample, pop in SIMULATED:
        haps, record = [], []
        for _ in range(2):
            hap, stats = simulate(chrom, morgans, src, rng)
            haps.append(hap)
            tot = [sum(v[i] for v in stats.values()) for i in range(3)]
            sub = [sum(stats[c][i] for c in (20, 21, 22)) for i in range(3)]
            record.append({"tracts": int(tot[0]), "french_fraction": tot[1] / tot[2],
                           "french_fraction_chr20_22": sub[1] / sub[2], "morgans": tot[2]})
        g = haps[0].astype(np.int16) + haps[1]
        g[(haps[0] == 9) | (haps[1] == 9)] = 9
        sim_geno.append(g.astype(np.uint8))
        samples.append({"id": sample, "pop": pop, "haplotypes": record,
                        "french_fraction": float(np.mean([h["french_fraction"] for h in record])),
                        "french_fraction_chr20_22": float(np.mean([h["french_fraction_chr20_22"] for h in record]))})

    order = [index[s] for p in ("Yoruba", "French") for s in reference[p]]
    out_inds = [(s, "U", p) for s, p in SIMULATED] + [(inds[k][0], inds[k][1], LABELS[inds[k][2]]) for k in order]
    out_geno = np.concatenate([np.stack(sim_geno, axis=1), geno[:, order]], axis=1)
    sub = np.flatnonzero(chrom >= 20)
    with tempfile.TemporaryDirectory() as tmp:
        write_packed(full, out_inds, snps, out_geno, Path(tmp))
        write_packed(small, out_inds, [snps[r] for r in sub], out_geno[sub], Path(tmp))

    truth = {
        "release": "AADR v66.p1, 1240K public (doi:10.7910/DVN/FFIDCW), CC0",
        "aadr": {ext: json.loads((CACHE / "manifest.json").read_text(encoding="utf-8"))["inputs"][ext]["md5"]
                 for ext in ("geno", "snp", "ind")},
        "snps": {"selection": f"every {THIN}rd autosomal SNP of the 1240K panel (rows 0, {THIN}, ...)",
                 "family_packed": len(snps), "family_chr20_22": len(sub)},
        "reference_samples": reference,
        "source_samples": sources,
        "duplicates_left_out": duplicates,
        # KING-robust, among the references and sources of each population
        "largest_kinship": largest,
        "simulation": {"seed": SEED, "numpy": np.__version__, "generations": GENERATIONS, "fraction": FRACTION,
                       "model": "crossovers at rate generations per Morgan on the AADR map; each segment French "
                                "with probability fraction, from a random source sample, one allele drawn per SNP"},
        "samples": samples,
        "outputs": {str(p.relative_to(out)) + f".{ext}": md5(Path(f"{p}.{ext}"))
                    for p in (small, full) for ext in ("geno", "snp", "ind")},
    }
    truth_path.write_text(json.dumps(truth, indent=1) + "\n", encoding="utf-8")
    for s in samples:
        print(f"{s['id']} {s['pop']}: French {s['french_fraction']:.3f} (chr20-22 {s['french_fraction_chr20_22']:.3f})")
    print(f"sources {sources}; duplicates left out {duplicates}; largest kinship {largest}")
    print(f"{len(out_inds)} samples; {len(snps)} SNPs, {len(sub)} on chr20-22")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else HERE)
