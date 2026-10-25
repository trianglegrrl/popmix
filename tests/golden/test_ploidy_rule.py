"""ADMIXTOOLS2-0007, live: a diploid genome whose first 1,000 SNPs hold no
heterozygous call.

admixtools2 calls a sample pseudo-haploid when it has no heterozygous call
among the first 1,000 SNPs. AADR's NA18606.DG (in CHB.DG), a high-coverage
genome, has 959 calls there and none heterozygous, so admixtools2 counts its
later heterozygous calls as half an allele. popmix's default looks at every
SNP; ``pseudohaploid="admixtools2"`` keeps the first-1,000 rule.

AADR v62.0 (popmix cannot read v66.p1's TGENO yet, task O1).
"""

from __future__ import annotations

from pathlib import Path

import pytest

import popmix as pm

pytestmark = pytest.mark.slow

AADR = "/references/aadr/v62.0_1240k_public"
POPS = ["CHB.DG", "Han.DG", "French.DG", "Mbuti.DG"]


@pytest.fixture(scope="module")
def rules():
    if not Path(f"{AADR}.geno").exists():
        pytest.skip("AADR v62.0 not available")
    return {rule: pm.extract_f2(AADR, POPS, pseudohaploid=rule) for rule in ("detect", "admixtools2", "diploid")}


def test_only_na18606_is_reclassified():
    if not Path(f"{AADR}.geno").exists():
        pytest.skip("AADR v62.0 not available")
    detect = pm.sample_ploidy(AADR, POPS)
    first1000 = pm.sample_ploidy(AADR, POPS, pseudohaploid="admixtools2")
    changed = {s for s in detect if detect[s] != first1000[s]}
    assert changed == {"NA18606.DG"}
    assert (detect["NA18606.DG"], first1000["NA18606.DG"]) == (2, 1)


def test_isolation_only_chb_statistics_change(rules):
    detect, compat = rules["detect"], rules["admixtools2"]
    for a, b in [("Han.DG", "French.DG"), ("Han.DG", "Mbuti.DG"), ("French.DG", "Mbuti.DG")]:
        assert detect.f2(a, b) == compat.f2(a, b), (a, b)
    for other in ("Han.DG", "French.DG", "Mbuti.DG"):
        assert detect.f2("CHB.DG", other)["estimate"] != compat.f2("CHB.DG", other)["estimate"], other


def test_oracle_the_default_treats_the_genome_as_diploid(rules):
    """Every sample of these populations has heterozygous calls, so the
    default must equal counting all of them as diploid."""
    detect, diploid = rules["detect"], rules["diploid"]
    for i, a in enumerate(POPS):
        for b in POPS[i + 1:]:
            assert detect.f2(a, b) == diploid.f2(a, b), (a, b)
