"""Without a genetic map, the 2 Mb fallback blocks are said (TASKS.md G5).

admixtools2 makes jackknife blocks of 2 Mb of base-pair position when every
SNP has the same map position; those are shorter than the 5 cM default, and
the standard errors come out slightly small.
"""

from __future__ import annotations

import shutil
import subprocess
import warnings
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = ROOT / "tests" / "golden" / "dates" / "family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]
MESSAGE = ("no genetic map (every SNP at the same map position): jackknife blocks are 2 Mb of base-pair position, "
           "as admixtools2, shorter than the 5 cM default, so standard errors come out slightly small; add genetic "
           "positions to the .snp or .bim file")


def caught(f) -> list[str]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f()
    return [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]


@pytest.fixture(scope="module")
def no_map(tmp_path_factory) -> str:
    prefix = tmp_path_factory.mktemp("nomap") / "data"
    for ext in ("geno", "ind"):
        shutil.copy(f"{SMALL}.{ext}", f"{prefix}.{ext}")
    rows = [line.split() for line in Path(f"{SMALL}.snp").read_text().splitlines()]
    Path(f"{prefix}.snp").write_text("".join(" ".join([r[0], r[1], "0", *r[3:]]) + "\n" for r in rows))
    return str(prefix)


def test_extract_f2_warns_without_a_map(no_map):
    assert MESSAGE in caught(lambda: pm.extract_f2(no_map, POPS))


def test_extract_f2_with_a_map_does_not_warn():
    assert MESSAGE not in caught(lambda: pm.extract_f2(str(SMALL), POPS))


def test_base_pair_blocks_asked_for_do_not_warn(no_map):
    assert MESSAGE not in caught(lambda: pm.extract_f2(no_map, POPS, block_size=2e6))


def test_genotype_statistics_warn_without_a_map(no_map):
    assert MESSAGE in caught(lambda: pm.qpdstat(no_map, "French", "Yoruba", "Simulation", "Simulation.1"))


def test_f2_from_frequencies_warns_without_a_map(no_map):
    af = pm.read_genotypes(no_map).allele_frequencies()
    assert MESSAGE in caught(lambda: pm.f2_from_frequencies(af))


def test_the_cli_says_so(no_map):
    out = subprocess.run(["popmix", "f2", "-i", no_map, "--pop1", "French", "--pop2", "Yoruba"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert f"Warning: {MESSAGE}" in out.stderr
