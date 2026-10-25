"""extract_f2 gives the same f2 blocks whatever format the genotypes are in.

PACKEDANCESTRYMAP input is streamed through the pipeline checked against
admixtools2; every other format must give the same blocks.
"""

import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def _convert(tmp_path: Path, fmt: str) -> str:
    out = tmp_path / f"small_{fmt}"
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", str(out), "-of", fmt], check=True,
                   capture_output=True)
    return str(out)


def _assert_same_blocks(ours, theirs):
    assert list(ours.block_lengths) == list(theirs.block_lengths)
    for i, a in enumerate(POPS):
        for b in POPS[i + 1:]:
            np.testing.assert_allclose(ours.f2_per_block(a, b), theirs.f2_per_block(a, b), rtol=1e-12, atol=0)


@pytest.mark.parametrize("fmt", ["packedancestrymap", "eigenstrat", "plink"])
def test_converted_genotypes_give_the_same_f2_blocks(tmp_path, fmt):
    original = pm.extract_f2(SMALL, POPS)
    converted = pm.extract_f2(_convert(tmp_path, fmt), POPS)
    _assert_same_blocks(converted, original)


def test_a_one_sample_population_matches_across_formats(tmp_path):
    """Simulation is a single sample, where the per-population bias correction matters most."""
    original = pm.extract_f2(SMALL, POPS)
    text = pm.extract_f2(_convert(tmp_path, "eigenstrat"), POPS)
    assert text.f2("Simulation", "French") == original.f2("Simulation", "French")


def test_statistic_dicts_have_a_fixed_key_order():
    """The result dicts (and so the CLI's JSON) do not depend on hash order."""
    f2 = pm.extract_f2(SMALL, POPS)
    results = [f2.f2("French", "Yoruba"), f2.f3("Simulation", "French", "Yoruba"),
               f2.f4("Simulation", "Simulation.1", "French", "Yoruba"), f2.fst("French", "Yoruba")]
    for r in results:
        assert list(r) == sorted(r)


def test_cli_json_is_byte_identical_across_runs():
    args = ["popmix", "f2", "--input", SMALL, "--pop1", "French", "--pop2", "Yoruba", "--format", "json"]
    outputs = {subprocess.run(args, check=True, capture_output=True, text=True).stdout for _ in range(4)}
    assert len(outputs) == 1


def _ind(prefix: str) -> list[tuple[str, str]]:
    return [(line.split()[0], line.split()[2]) for line in Path(f"{prefix}.ind").read_text().splitlines()]


@pytest.mark.parametrize("fmt", ["packedancestrymap", "eigenstrat"])
def test_inds_and_pops_select_and_relabel_samples(tmp_path, fmt):
    """As admixtools2 f2_from_geno(inds, pops): only the listed samples, with these labels."""
    prefix = SMALL if fmt == "packedancestrymap" else _convert(tmp_path, fmt)
    rows = [r for r in _ind(SMALL) if r[1] in POPS]
    inds, pops = [r[0] for r in rows], [r[1] for r in rows]
    _assert_same_blocks(pm.extract_f2(prefix, POPS, inds=inds, pops=pops), pm.extract_f2(SMALL, POPS))

    # leaving out a French sample changes f2 with French, and nothing else
    rest = [k for k, p in enumerate(pops) if p == "French"][0]
    keep = [k for k in range(len(inds)) if k != rest]
    fewer = pm.extract_f2(prefix, POPS, inds=[inds[k] for k in keep], pops=[pops[k] for k in keep])
    assert fewer.f2("French", "Yoruba") != pm.extract_f2(SMALL, POPS).f2("French", "Yoruba")

    # relabelled samples form a new population
    relabelled = ["Half" if p == "French" and k % 2 else p for k, p in enumerate(pops)]
    f2 = pm.extract_f2(prefix, POPS + ["Half"], inds=inds, pops=relabelled)
    assert sorted(f2.populations) == sorted(POPS + ["Half"])
    assert f2.f2("Half", "French") is not None


def test_inds_and_pops_must_match():
    with pytest.raises(ValueError, match="same length"):
        pm.extract_f2(SMALL, POPS, inds=["a", "b"], pops=["French"])
    with pytest.raises(ValueError, match="twice"):
        pm.extract_f2(SMALL, POPS, inds=["a", "a"], pops=["French", "French"])
    with pytest.raises(ValueError, match="not in"):
        pm.extract_f2(SMALL, ["French"], inds=["nobody"], pops=["French"])
