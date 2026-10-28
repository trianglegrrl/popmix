"""SNP filters for f2 extraction, in Python and on the CLI (TASKS.md H1).

The SNP counts are admixtools2 2.0.10's f2_from_geno on the same data and
populations (polymorphic SNPs in its blocks); the AADR goldens are
tests/golden/test_snp_filters_admixtools2.py.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]
IDS = [line.split()[0] for line in Path(f"{SMALL}.snp").read_text().splitlines()]


def n_snps(**kw) -> int:
    return sum(pm.extract_f2(SMALL, POPS, **kw).block_lengths)


@pytest.mark.parametrize(("kw", "admixtools2"), [
    ({}, 21080),
    ({"transitions": False}, 254),
    ({"transversions": False}, 20826),
    ({"minmaf": 0.1}, 17647),
    ({"maxmaf": 0.2}, 8666),
    ({"minmaf": 0.05, "maxmaf": 0.3}, 11998),
])
def test_snp_counts_are_admixtools2s(kw, admixtools2):
    assert n_snps(**kw) == admixtools2


def test_keepsnps_overrides_the_other_filters():
    keep = pm.extract_f2(SMALL, POPS, keepsnps=IDS[:5000])
    assert pm.extract_f2(SMALL, POPS, keepsnps=IDS[:5000], minmaf=0.1, transitions=False).block_lengths == \
        keep.block_lengths


def test_excluding_all_but_some_snps_keeps_those():
    keep = pm.extract_f2(SMALL, POPS, keepsnps=IDS[:5000])
    rest = pm.extract_f2(SMALL, POPS, excludesnps=IDS[5000:])
    assert rest.block_lengths == keep.block_lengths
    assert rest.f2("French", "Yoruba") == keep.f2("French", "Yoruba")


def test_a_maf_bound_past_a_half_is_an_error():
    with pytest.raises(pm.InputError, match="minmaf"):
        pm.extract_f2(SMALL, POPS, minmaf=0.6)


def popmix(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["popmix", *args], capture_output=True, text=True)


def test_extract_f2_takes_transversions_only(tmp_path):
    r = popmix("extract-f2", "-i", SMALL, "--populations", ",".join(POPS), "--transversions-only",
               "-o", str(tmp_path / "f2"))
    assert r.returncode == 0, r.stderr
    f2 = pm.read_f2(str(tmp_path / "f2"))
    assert sum(f2.block_lengths) == 254
    assert f2.provenance["parameters"]["transitions"] is False
    sidecar = json.loads((tmp_path / "f2" / "provenance.json").read_text())
    assert sidecar["parameters"]["transversions_only"] is True


def test_statistics_on_genotypes_take_the_filters(tmp_path):
    keep = tmp_path / "keep.snp"
    keep.write_text("".join(line + "\n" for line in Path(f"{SMALL}.snp").read_text().splitlines()[:5000]))
    r = popmix("f2", "--pop1", "French", "--pop2", "Yoruba", "-i", SMALL, "--keep-snps", str(keep),
               "--minmaf", "0.05", "--format", "json")
    assert r.returncode == 0, r.stderr
    want = pm.extract_f2(SMALL, ["French", "Yoruba"], keepsnps=IDS[:5000], minmaf=0.05).f2("French", "Yoruba")
    assert json.loads(r.stdout)["est"] == pytest.approx(want["est"], rel=1e-12)


def test_qpadm_takes_the_filters():
    r = popmix("qpadm", "-i", SMALL, "--target", "Simulation", "--left", "French", "--right",
               "Yoruba,Simulation.1", "--maxmaf", "0.2", "--format", "json")
    assert r.returncode == 0, r.stderr
    f2 = pm.extract_f2(SMALL, ["Simulation", "French", "Yoruba", "Simulation.1"], maxmaf=0.2)
    want = pm.qpadm(f2, ["French"], ["Yoruba", "Simulation.1"], "Simulation")
    assert json.loads(r.stdout)["p"] == pytest.approx(want.p, rel=1e-10)


def test_both_only_flags_are_a_usage_error():
    r = popmix("extract-f2", "-i", SMALL, "--transversions-only", "--transitions-only", "-o", "x")
    assert r.returncode == 2
    assert "leave no SNPs" in r.stderr
