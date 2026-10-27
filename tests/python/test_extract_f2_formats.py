"""extract_f2 gives the same f2 blocks whatever format the genotypes are in.

PACKEDANCESTRYMAP input is streamed through the pipeline checked against
admixtools2; every other format must give the same blocks.
"""

import json
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
    # f2 block sums are matrix products: the formats agree to rounding
    got, want = text.f2("Simulation", "French"), original.f2("Simulation", "French")
    assert got.keys() == want.keys()
    for key in want:
        assert got[key] == pytest.approx(want[key], rel=1e-12), key


def test_statistic_dicts_have_a_fixed_key_order():
    """The result dicts (and so the CLI's JSON) have the keys of a table row,
    in its order (TASKS.md J2), whatever the hash order."""
    f2 = pm.extract_f2(SMALL, POPS)
    assert list(f2.f2("French", "Yoruba")) == ["est", "se", "z", "p", "n_snps"]
    assert list(f2.f3("Simulation", "French", "Yoruba")) == ["est", "se", "z", "p", "n_blocks"]
    assert list(f2.f4("Simulation", "Simulation.1", "French", "Yoruba")) == ["est", "se", "z", "p", "n_blocks"]
    assert list(f2.fst("French", "Yoruba")) == ["est", "se", "n_snps"]
    assert isinstance(f2.f2("French", "Yoruba")["n_snps"], int)


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


@pytest.fixture(scope="module")
def pseudohaploid(tmp_path_factory) -> str:
    """The small data with Simulation.1 made pseudohaploid: each het call
    becomes a homozygous one, as in ancient-DNA pseudohaploid calls."""
    d = tmp_path_factory.mktemp("pseudohaploid")
    text = d / "text"
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", str(text), "-of", "eigenstrat"], check=True,
                   capture_output=True)
    inds = Path(f"{text}.ind").read_text().splitlines()
    cols = [i for i, line in enumerate(inds) if line.split()[2] == "Simulation.1"]
    rows = []
    for j, line in enumerate(Path(f"{text}.geno").read_text().splitlines()):
        row = list(line)
        for i in cols:
            if row[i] == "1":
                row[i] = "2" if (i + j) % 2 else "0"
        rows.append("".join(row))
    Path(f"{text}.geno").write_text("\n".join(rows) + "\n")
    return str(text)


def _from(prefix: str, tmp_path: Path, fmt: str) -> str:
    if fmt == "eigenstrat":
        return prefix
    out = tmp_path / fmt
    subprocess.run(["popmix", "convertf", "-i", prefix, "-o", str(out), "-of", fmt], check=True,
                   capture_output=True)
    return str(out)


@pytest.mark.parametrize("fmt", ["eigenstrat", "plink", "ped"])
def test_pseudohaploid_samples_are_detected_in_every_format(pseudohaploid, tmp_path, fmt):
    """Packed input detects pseudohaploid samples (no het in the first 1,000
    SNPs) as admixtools2 does; every other format must too."""
    packed = pm.extract_f2(_from(pseudohaploid, tmp_path, "packedancestrymap"), POPS)
    ours = pm.extract_f2(_from(pseudohaploid, tmp_path, fmt), POPS)
    _assert_same_blocks(ours, packed)


def test_pseudohaploid_samples_are_detected_when_reading_genotypes(pseudohaploid, tmp_path):
    """README "Method 1": read_genotypes, allele_frequencies, f2_from_frequencies."""
    packed = _from(pseudohaploid, tmp_path, "packedancestrymap")
    want = pm.extract_f2(packed, POPS).f2("Simulation.1", "French")["est"]
    for prefix in (packed, pseudohaploid):
        af = pm.read_genotypes(prefix).allele_frequencies(populations=POPS)
        got = pm.f2_from_frequencies(af).f2("Simulation.1", "French")["est"]
        np.testing.assert_allclose(got, want, rtol=1e-12, err_msg=prefix)


@pytest.fixture(scope="module")
def late_hets(tmp_path_factory) -> str:
    """The small data with French's first sample homozygous over the first
    1,000 SNPs: diploid (hets later on), but admixtools2's rule (no het in
    the first 1,000 SNPs) calls it pseudohaploid."""
    d = tmp_path_factory.mktemp("late_hets")
    text = d / "text"
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", str(text), "-of", "eigenstrat"], check=True,
                   capture_output=True)
    inds = Path(f"{text}.ind").read_text().splitlines()
    col = next(i for i, line in enumerate(inds) if line.split()[2] == "French")
    rows = Path(f"{text}.geno").read_text().splitlines()
    for j in range(1000):
        if rows[j][col] == "1":
            rows[j] = rows[j][:col] + "2" + rows[j][col + 1:]
    assert any(row[col] == "1" for row in rows[1000:])
    Path(f"{text}.geno").write_text("\n".join(rows) + "\n")
    return str(text)


@pytest.mark.parametrize("fmt", ["eigenstrat", "packedancestrymap"])
def test_a_sample_with_any_het_call_is_diploid(late_hets, tmp_path, fmt):
    """ADMIXTOOLS2-0007: by default a sample is pseudohaploid only with no het
    call at all; pseudohaploid="admixtools2" keeps the first-1,000-SNP rule."""
    prefix = _from(late_hets, tmp_path, fmt)
    default = pm.extract_f2(prefix, POPS)
    diploid = pm.extract_f2(prefix, POPS, pseudohaploid="diploid")
    compat = pm.extract_f2(prefix, POPS, pseudohaploid="admixtools2")
    _assert_same_blocks(default, diploid)
    assert compat.f2("French", "Yoruba")["est"] != default.f2("French", "Yoruba")["est"]
    # isolation: statistics without French are unaffected, once the SNP set
    # does not depend on French's frequencies (poly_only keeps SNPs whose
    # frequencies differ, and recounting French moves them)
    fixed = {rule: pm.extract_f2(prefix, POPS, poly_only=False, pseudohaploid=rule)
             for rule in ("detect", "admixtools2")}
    assert fixed["detect"].f2("Yoruba", "Simulation") == fixed["admixtools2"].f2("Yoruba", "Simulation")
    assert fixed["detect"].f2("French", "Yoruba") != fixed["admixtools2"].f2("French", "Yoruba")


def test_pseudohaploid_rule_names_are_checked():
    with pytest.raises(ValueError, match="pseudohaploid"):
        pm.extract_f2(SMALL, POPS, pseudohaploid="first1000")


def test_cli_pseudohaploid_option(late_hets):
    def f2(*extra):
        out = subprocess.run(["popmix", "f2", "--input", late_hets, "--pop1", "French", "--pop2", "Yoruba",
                              "--format", "json", *extra], check=True, capture_output=True, text=True)
        return json.loads(out.stdout)["est"]

    assert f2() == f2("--pseudohaploid", "diploid")
    assert f2() != f2("--pseudohaploid", "admixtools2")


def test_text_with_crlf_lines_gives_the_same_f2_blocks(tmp_path):
    """Text .geno of even lines is read in place; other text goes through the reader."""
    text = _convert(tmp_path, "eigenstrat")
    geno = Path(f"{text}.geno")
    geno.write_bytes(geno.read_bytes().replace(b"\n", b"\r\n"))
    _assert_same_blocks(pm.extract_f2(text, POPS), pm.extract_f2(SMALL, POPS))
