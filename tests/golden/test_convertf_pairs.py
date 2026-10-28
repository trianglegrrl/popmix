"""popmix convertf against AdmixTools 8 convertf on every format pair (TASKS.md P5, K3).

The small family data is written by convertf in each of its six formats;
each of those is then converted to each format by convertf and by
``popmix convertf``, and the two outputs must hold the same samples, SNPs
and genotypes. PED and PACKEDPED carry no population labels (convertf
writes a family number), so from those inputs and to those outputs the
samples are compared by ID; convertf names a sample read from them
FID:IID, popmix IID (DEVIATIONS.md). Genotypes are compared as allele
letters: a PED .map without alleles (convertf writes one) leaves each tool
to choose the counted allele.
"""

from __future__ import annotations

import itertools
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

SMALL = Path(__file__).parent / "dates" / "family_chr20_22"
CONVERTF = os.environ.get("POPMIX_CONVERTF",
                          str(Path.home() / ".local/share/popmix-references/admixtools-8.0.2/bin/convertf"))
pytestmark = [pytest.mark.slow,
              pytest.mark.skipif(not Path(CONVERTF).exists(), reason="AdmixTools 8 convertf not installed")]

# convertf name -> (popmix name, the extensions popmix reads and writes)
FORMATS = {
    "PACKEDANCESTRYMAP": ("packedancestrymap", ("geno", "snp", "ind")),
    "EIGENSTRAT": ("eigenstrat", ("geno", "snp", "ind")),
    "ANCESTRYMAP": ("ancestrymap", ("ancestrymapgeno", "snp", "ind")),
    "PACKEDPED": ("plink", ("bed", "bim", "fam")),
    "PED": ("ped", ("ped", "map", "pedind")),
    "TRANSPOSE_PACKED": ("tgeno", ("geno", "snp", "ind")),
}
UNLABELLED = {"PACKEDPED", "PED"}


def convertf(src: Path, src_fmt: str, out: Path, out_fmt: str) -> None:
    (g, s, i), (og, os_, oi) = FORMATS[src_fmt][1], FORMATS[out_fmt][1]
    par = out.with_suffix(".par")
    par.write_text(f"genotypename: {src}.{g}\nsnpname: {src}.{s}\nindivname: {src}.{i}\noutputformat: {out_fmt}\n"
                   f"genotypeoutname: {out}.{og}\nsnpoutname: {out}.{os_}\nindivoutname: {out}.{oi}\n")
    subprocess.run([CONVERTF, "-p", str(par)], check=True, capture_output=True, cwd=out.parent)


@pytest.fixture(scope="module")
def inputs(tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("inputs")
    for fmt in FORMATS:
        convertf(SMALL, "PACKEDANCESTRYMAP", d / fmt, fmt)
    return d


def _letters(g) -> np.ndarray:
    """Each genotype as its two alleles, in sorted order ("" when missing)."""
    codes = g.to_numpy()
    first = np.array([s.allele_ref for s in g.snps])[:, None]
    second = np.array([s.allele_alt for s in g.snps])[:, None]
    swap = first > second
    lo, hi = np.where(swap, second, first), np.where(swap, first, second)
    het = np.char.add(lo, hi)
    return np.select([codes == 2, codes == 1, codes == 0],
                     [np.char.add(first, first), np.broadcast_to(het, codes.shape), np.char.add(second, second)], "")


def _iid(sample_id: str) -> str:
    return sample_id.split(":", 1)[1]


def _same(ours: Path, theirs: Path, src_fmt: str, out_fmt: str) -> None:
    a, b = pm.read_genotypes(str(ours)), pm.read_genotypes(str(theirs))
    if src_fmt in UNLABELLED:
        assert [s.id for s in a.samples] == [_iid(s.id) for s in b.samples]
    elif out_fmt in UNLABELLED:
        assert [s.id for s in a.samples] == [s.id for s in b.samples]
    else:
        assert [(s.id, s.population) for s in a.samples] == [(s.id, s.population) for s in b.samples]
    assert [(s.id, s.chromosome, s.position_bp) for s in a.snps] == \
        [(s.id, s.chromosome, s.position_bp) for s in b.snps]
    np.testing.assert_array_equal(_letters(a), _letters(b))


@pytest.mark.parametrize(("src_fmt", "out_fmt"), list(itertools.product(FORMATS, FORMATS)))
def test_popmix_writes_what_convertf_writes(inputs, tmp_path, src_fmt, out_fmt):
    theirs, ours = tmp_path / "theirs", tmp_path / "ours"
    convertf(inputs / src_fmt, src_fmt, theirs, out_fmt)
    r = subprocess.run(["popmix", "convertf", "-i", str(inputs / src_fmt), "-if", FORMATS[src_fmt][0],
                        "-o", str(ours), "-of", FORMATS[out_fmt][0]], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    _same(ours, theirs, src_fmt, out_fmt)
