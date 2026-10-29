"""convertf parameter files against AdmixTools 8 convertf (TASKS.md P2).

Each parameter file runs through convertf and through ``popmix convertf
-p``; the two outputs must hold the same samples, SNPs and genotypes. The
files use a DIR name, input files of different stems, poplistname,
badsnpname (with a SNP the data does not have), chrom with lopos/hipos,
both spellings of the output keys, and the default output format.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

from .genotype_letters import letters

ROOT = Path(__file__).parents[1].parent
SMALL = ROOT / "tests" / "golden" / "dates" / "family_chr20_22"
CONVERTF = os.environ.get("POPMIX_CONVERTF",
                          str(Path.home() / ".local/share/popmix-references/admixtools-8.0.2/bin/convertf"))
pytestmark = pytest.mark.skipif(not Path(CONVERTF).exists(), reason="AdmixTools 8 convertf not installed")

OUT = {"PACKEDPED": ("bed", "bim", "fam"), "PED": ("ped", "map", "pedind")}


def write_par(path: Path, out: Path, fmt: str | None, extra: dict[str, str], mergeit_keys: bool = False) -> Path:
    geno, snp, ind = OUT.get(fmt or "", ("geno", "snp", "ind"))
    keys = ("genooutfilename", "snpoutfilename", "indoutfilename") if mergeit_keys else \
        ("genotypeoutname", "snpoutname", "indivoutname")
    lines = {"DIR": str(SMALL.parent), "genotypename": f"DIR/{SMALL.name}.geno",
             "snpname": f"DIR/{SMALL.name}.snp", "indivname": str(path.parent / "renamed.ind"),
             **({"outputformat": fmt} if fmt else {}),
             keys[0]: f"{out}.{geno}", keys[1]: f"{out}.{snp}", keys[2]: f"{out}.{ind}", **extra}
    path.write_text("".join(f"{k}: {v}\n" for k, v in lines.items()), encoding="utf-8")
    return path


@pytest.fixture
def inputs(tmp_path) -> Path:
    (tmp_path / "renamed.ind").write_text(Path(f"{SMALL}.ind").read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "pops.txt").write_text("Yoruba\nSimulation.1\n", encoding="utf-8")
    snps = [line.split()[0] for line in Path(f"{SMALL}.snp").read_text(encoding="utf-8").splitlines()]
    (tmp_path / "bad.snp").write_text("".join(f"{s}\n" for s in snps[::7]) + "rs_not_in_the_data\n", encoding="utf-8")
    return tmp_path


def run_both(inputs: Path, fmt: str | None, extra: dict[str, str], mergeit_keys: bool = False):
    ref = write_par(inputs / "ref.par", inputs / "ref", fmt, extra, mergeit_keys)
    ours = write_par(inputs / "ours.par", inputs / "ours", fmt, extra, mergeit_keys)
    subprocess.run([CONVERTF, "-p", str(ref)], check=True, capture_output=True, cwd=inputs)
    r = subprocess.run(["popmix", "convertf", "-p", str(ours)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return inputs


def same(d: Path, fmt: str | None) -> None:
    """popmix's output (d/ours) holds what convertf's (d/ref) does."""
    if fmt == "ANCESTRYMAP":
        # SNP, sample, genotype lines; read as text, both named .geno
        def lines(name: str) -> list[list[str]]:
            return [line.split() for line in Path(d / name).read_text(encoding="utf-8").splitlines()]
        for ext in ("geno", "snp", "ind"):
            assert lines(f"ours.{ext}") == lines(f"ref.{ext}"), ext
        return
    a, b = pm.read_genotypes(str(d / "ours")), pm.read_genotypes(str(d / "ref"))
    # convertf's .fam and .pedind have a family number where popmix writes the population
    key = (lambda s: s.id) if fmt in ("PACKEDPED", "PED") else (lambda s: (s.id, s.population))
    assert [key(s) for s in a.samples] == [key(s) for s in b.samples]
    assert [(s.id, s.chromosome, s.position_bp) for s in a.snps] == [(s.id, s.chromosome, s.position_bp)
                                                                      for s in b.snps]
    if fmt == "PED":
        # convertf's .map has no alleles
        np.testing.assert_array_equal(letters(a), letters(b))
    else:
        np.testing.assert_array_equal(a.to_numpy(), b.to_numpy())


@pytest.mark.parametrize("fmt", [None, "EIGENSTRAT", "PACKEDPED", "PED", "ANCESTRYMAP", "TRANSPOSE_PACKED"])
def test_populations_and_bad_snps(inputs, fmt):
    d = run_both(inputs, fmt, {"poplistname": str(inputs / "pops.txt"), "badsnpname": str(inputs / "bad.snp")})
    same(d, fmt)


def test_a_region_with_mergeit_key_names(inputs):
    d = run_both(inputs, "EIGENSTRAT", {"chrom": "21", "lopos": "20000000", "hipos": "35000000"},
                 mergeit_keys=True)
    same(d, "EIGENSTRAT")
    assert {s.chromosome for s in pm.read_genotypes(str(d / "ref")).snps} == {"21"}


def test_an_unsupported_key_is_an_error(inputs):
    par = write_par(inputs / "x.par", inputs / "x", "EIGENSTRAT", {"killr2": "YES"})
    r = subprocess.run(["popmix", "convertf", "-p", str(par)], capture_output=True, text=True)
    assert r.returncode == 1
    assert "does not support convertf's killr2" in r.stderr


def test_popmix_reads_convertfs_ancestrymap(inputs):
    """ANCESTRYMAP counts copies of the first allele and writes -1 for a missing call; popmix
    read and wrote it inverted, and refused -1, until this test."""
    d = run_both(inputs, "ANCESTRYMAP", {})
    (d / "ref.geno").rename(d / "ref.ancestrymapgeno")
    out = subprocess.run(["popmix", "convertf", "-i", str(d / "ref"), "-if", "ancestrymap", "-o", str(d / "back"),
                          "-of", "eigenstrat"], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    np.testing.assert_array_equal(pm.read_genotypes(str(d / "back")).to_numpy(), pm.read_genotypes(str(SMALL)).to_numpy())
