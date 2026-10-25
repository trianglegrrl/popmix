"""ADMIXTOOLS2-0007, live: diploid genomes whose first 1,000 SNPs hold no
heterozygous call.

admixtools2 calls a sample pseudo-haploid when it has no heterozygous call
among the first 1,000 SNPs. AADR's NA18606.DG (in CHB.DG) and HGDP00581.DG
(in Druze.DG, and so in the dg_11pop goldens), high-coverage genomes, have
no heterozygous call there, so admixtools2 counts their later heterozygous
calls as half an allele. popmix's default looks at every SNP;
``pseudohaploid="admixtools2"`` keeps the first-1,000 rule.

AADR v66.p1, labelled by data type (tests/golden/aadr/build_aadr_v66p1.py).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import popmix as pm

pytestmark = pytest.mark.slow

AADR = str(Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references"))
           / "aadr_v66p1" / "aadr_v66p1")
OTHERS = ["Han.DG", "French.DG", "Mbuti.DG"]
# the genome, and its population
CASES = {"NA18606.DG": "CHB.DG", "HGDP00581.DG": "Druze.DG"}
_rules: dict = {}


def _pops(sample: str) -> list[str]:
    return [CASES[sample], *OTHERS]


def rules(sample: str) -> dict:
    if not Path(f"{AADR}.geno").exists():
        pytest.skip("AADR v66.p1 dataset not built")
    if sample not in _rules:
        _rules[sample] = {rule: pm.extract_f2(AADR, _pops(sample), pseudohaploid=rule)
                          for rule in ("detect", "admixtools2", "diploid")}
    return _rules[sample]


@pytest.mark.parametrize("sample", CASES)
def test_only_that_genome_is_reclassified(sample):
    if not Path(f"{AADR}.geno").exists():
        pytest.skip("AADR v66.p1 dataset not built")
    detect = pm.sample_ploidy(AADR, _pops(sample))
    first1000 = pm.sample_ploidy(AADR, _pops(sample), pseudohaploid="admixtools2")
    changed = {s for s in detect if detect[s] != first1000[s]}
    assert changed == {sample}
    assert (detect[sample], first1000[sample]) == (2, 1)


@pytest.mark.parametrize("sample", CASES)
def test_isolation_only_its_population_changes(sample):
    r = rules(sample)
    detect, compat = r["detect"], r["admixtools2"]
    for i, a in enumerate(OTHERS):
        for b in OTHERS[i + 1:]:
            assert detect.f2(a, b) == compat.f2(a, b), (a, b)
    for other in OTHERS:
        assert detect.f2(CASES[sample], other)["estimate"] != compat.f2(CASES[sample], other)["estimate"], other


@pytest.mark.parametrize("sample", CASES)
def test_oracle_the_default_treats_the_genome_as_diploid(sample):
    """Every sample of these populations has heterozygous calls, so the
    default must equal counting all of them as diploid."""
    r = rules(sample)
    detect, diploid = r["detect"], r["diploid"]
    pops = _pops(sample)
    for i, a in enumerate(pops):
        for b in pops[i + 1:]:
            assert detect.f2(a, b) == diploid.f2(a, b), (a, b)
