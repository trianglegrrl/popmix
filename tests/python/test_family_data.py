"""The small family data: AADR v66.p1 genomes and admixed samples popmix simulates.

tests/golden/dates/build_family.py makes tests/golden/dates/family_chr20_22
(in the repository) and example_data/family_packed (every autosome, built
locally) from AADR v66.p1 (CC0), and records the simulation and its truth
in tests/golden/dates/family_truth.json. The reproducibility check, which
reruns the script on AADR, is tests/golden/test_family_build.py.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

ROOT = Path(__file__).resolve().parents[2]
DATES = ROOT / "tests" / "golden" / "dates"
SMALL = DATES / "family_chr20_22"
FULL = DATES / "example_data" / "family_packed"
TRUTH = json.loads((DATES / "family_truth.json").read_text(encoding="utf-8"))
SIMULATED = ["Simulation", "Simulation.1"]


def ind(prefix: Path) -> list[tuple[str, str, str]]:
    return [tuple(line.split()) for line in Path(f"{prefix}.ind").read_text(encoding="utf-8").splitlines()]


def test_composition():
    rows = ind(SMALL)
    assert Counter(pop for _, _, pop in rows) == {"French": 24, "Yoruba": 21, "Simulation": 1, "Simulation.1": 4}
    assert [s for s, _, _ in rows[:5]] == ["NA00000", "NA00001", "NA00002", "NA00003", "NA00004"]


def test_real_samples_are_the_recorded_aadr_genomes():
    rows = ind(SMALL)
    for pop in ("French", "Yoruba"):
        assert [s for s, _, p in rows if p == pop] == TRUTH["reference_samples"][pop]
        assert all(s.endswith(".DG") for s in TRUTH["reference_samples"][pop])


def test_sources_are_not_reference_samples():
    for pop in ("French", "Yoruba"):
        sources = set(TRUTH["source_samples"][pop])
        assert sources
        assert sources.isdisjoint(TRUTH["reference_samples"][pop])
    # SGDP genomes of HGDP individuals: one person, sequenced again
    again = {pair.split(" = ")[0] for pair in TRUTH["duplicates_left_out"]}
    assert len(again) == 6
    assert again.isdisjoint(sample for sample, _, _ in ind(SMALL))


def test_source_pools_and_relatedness():
    assert len(TRUTH["source_samples"]["French"]) == 4
    # AADR YRI.DG (1000 Genomes Yoruba, and HGDP00927), unrelated to each other and to the HGDP Yoruba
    assert len(TRUTH["source_samples"]["Yoruba"]) == 30
    assert all(s.endswith(".DG") for s in TRUTH["source_samples"]["Yoruba"])
    # no third-degree relatives or closer (KING-robust 2^-4.5)
    assert TRUTH["largest_kinship"]["French"] < 2 ** -4.5
    assert TRUTH["largest_kinship"]["Yoruba"] < 2 ** -4.5


def test_simulated_samples_test_as_admixed():
    """f3(Simulation.1; French, Yoruba) is negative: the sources' own drift does
    not hide the admixture."""
    f3 = pm.qp3pop(str(SMALL), ["Simulation.1"], ["French"], ["Yoruba"])
    assert f3["z"].iloc[0] < -3


def test_simulation_parameters():
    sim = TRUTH["simulation"]
    assert sim["generations"] == 50
    assert sim["fraction"] == {"French": 0.2, "Yoruba": 0.8}
    assert isinstance(sim["seed"], int)


def test_chromosomes_20_to_22_only():
    chroms = {line.split()[1] for line in Path(f"{SMALL}.snp").read_text(encoding="utf-8").splitlines()}
    assert chroms == {"20", "21", "22"}
    assert len(Path(f"{SMALL}.snp").read_text(encoding="utf-8").splitlines()) == TRUTH["snps"]["family_chr20_22"]


def test_realized_ancestry_is_near_the_simulated_fraction():
    for sample in TRUTH["samples"]:
        for hap in sample["haplotypes"]:
            assert 0.0 <= hap["french_fraction"] <= 1.0
            assert hap["tracts"] > 22
    mean = np.mean([h["french_fraction"] for s in TRUTH["samples"] for h in s["haplotypes"]])
    # 10 haplotypes of ~36 Morgans in tracts of 2 cM: the mean's SD is about 0.01
    assert mean == pytest.approx(0.2, abs=0.04)


def test_qpadm_recovers_the_truth_on_chr20_22(tmp_path):
    """qpAdm (f4-based, so unbiased by the reference samples' finite size)
    finds Simulation.1's realized French ancestry on chromosomes 20-22, with
    French and Yoruba each split in three by sample order (one of each as a
    source, the others as right populations)."""
    seen: Counter = Counter()
    rows = []
    for sample, sex, pop in ind(SMALL):
        if pop in ("French", "Yoruba"):
            seen[pop] += 1
            pop = f"{pop}_{seen[pop] % 3 + 1}"
        rows.append(f"{sample} {sex} {pop}\n")
    prefix = tmp_path / "family6"
    for ext in ("geno", "snp"):
        (tmp_path / f"family6.{ext}").write_bytes(Path(f"{SMALL}.{ext}").read_bytes())
    (tmp_path / "family6.ind").write_text("".join(rows), encoding="utf-8")
    fit = pm.qpadm(str(prefix), ["French_1", "Yoruba_1"], ["French_2", "Yoruba_2", "French_3", "Yoruba_3"],
                   "Simulation.1")
    truth = np.mean([s["french_fraction_chr20_22"] for s in TRUTH["samples"] if s["pop"] == "Simulation.1"])
    assert abs(fit.weights[0] - truth) < 3 * fit.weight_se[0]


def test_simulated_samples_are_diploid():
    ploidy = pm.sample_ploidy(str(SMALL), SIMULATED)
    assert set(ploidy.values()) == {2}


@pytest.mark.skipif(not Path(f"{FULL}.geno").exists(), reason="example_data not built (build_family.py)")
def test_chr20_22_is_the_full_genome_restricted():
    full = Path(f"{FULL}.snp").read_text(encoding="utf-8").splitlines()
    assert [line for line in full if line.split()[1] in ("20", "21", "22")] == \
        Path(f"{SMALL}.snp").read_text(encoding="utf-8").splitlines()
    assert ind(FULL) == ind(SMALL)
    assert sum(1 for _ in full) == TRUTH["snps"]["family_packed"]
    assert {line.split()[1] for line in full} == {str(c) for c in range(1, 23)}
