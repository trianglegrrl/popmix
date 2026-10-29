"""DATES says when a source has few observed alleles per SNP (TASKS.md G10).

DATES uses a SNP once each source has one called sample, so a source of a
few samples, or of pseudo-haploid samples with much missing data, gives
noisy weights; popmix reports each source's mean observed alleles per SNP
used and warns below 10.
"""

from __future__ import annotations

import json
import subprocess
import warnings
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")


def caught(f) -> list[str]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f()
    return [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]


def test_a_source_of_one_sample_warns():
    messages = caught(lambda: pm.dates(SMALL, "Simulation.1", ["French", "Simulation"], jackknife=False))
    assert messages == ["DATES: source Simulation has 2.0 observed alleles per SNP used (fewer than 10; a SNP needs "
                        "only one called sample): its frequencies, the SNP weights, are noisy, which weakens the "
                        "signal and widens the SE"]


def test_sources_with_many_samples_do_not_warn():
    assert caught(lambda: pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"], jackknife=False)) == []


def test_the_result_holds_each_sources_mean_alleles():
    r = pm.dates(SMALL, "Simulation.1", ["French", "Simulation"], jackknife=False)
    assert r.source_alleles[1] == 2.0
    # 24 French genomes, a few calls missing
    assert r.source_alleles[0] == pytest.approx(47.767, abs=1e-3)


def test_the_cli_json_has_them():
    out = subprocess.run(["popmix", "dates", "-i", SMALL, "--admixed", "Simulation.1", "--sources", "French,Simulation",
                          "--no-jackknife", "--format", "json"], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert json.loads(out.stdout)["source_alleles"][1] == 2.0
    assert "Warning: DATES: source Simulation has 2.0 observed alleles per SNP used" in out.stderr
