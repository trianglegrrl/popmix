"""f2 directories and result files record how they were made (TASKS.md H3).

An f2 directory's ``metadata.json`` carries the popmix version, the
extraction's parameters, the SNP counts and each sample's ploidy; every
command that writes files leaves a ``.provenance.json`` beside them with the
command, every parameter (defaults included) and the inputs.
"""

from __future__ import annotations

import json
import pickle
import subprocess
from pathlib import Path

import pytest

import popmix as pm
from popmix.cli import cli

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")
POPS = ["French", "Yoruba", "Simulation.1"]


@pytest.fixture(scope="module")
def f2() -> pm.F2Blocks:
    return pm.extract_f2(SMALL, POPS, block_size=0.02, maxmiss=0.1)


def samples_of(pops: list[str]) -> list[tuple[str, str]]:
    rows = (line.split() for line in Path(f"{SMALL}.ind").read_text().splitlines())
    return [(i, pop) for i, _, pop in rows if pop in pops]


def test_extract_f2_records_its_parameters(f2):
    p = f2.provenance
    assert p["popmix_version"] == pm.__version__
    assert p["source"] == SMALL
    assert p["parameters"] == {"populations": POPS, "auto_only": True, "block_size": 0.02, "maxmiss": 0.1,
                               "poly_only": True, "pseudohaploid": "detect", "relabelled_samples": 0}


def test_extract_f2_records_the_snps_used(f2):
    assert f2.provenance["n_snps_used"] == sum(f2.block_lengths)
    n_snps = len(Path(f"{SMALL}.snp").read_text().splitlines())
    assert f2.provenance["n_snps_read"] == n_snps


def test_extract_f2_records_each_samples_ploidy(f2):
    samples = f2.provenance["samples"]
    assert [(s["id"], s["population"]) for s in samples] == samples_of(POPS)
    assert {s["ploidy"] for s in samples if s["population"] != "Simulation.1"} == {"diploid"}


def test_the_ploidy_rule_shows_in_the_ploidy():
    f2 = pm.extract_f2(SMALL, POPS, pseudohaploid="diploid")
    assert {s["ploidy"] for s in f2.provenance["samples"]} == {"diploid"}
    assert f2.provenance["parameters"]["pseudohaploid"] == "diploid"


def test_provenance_survives_saving_reading_and_pickling(f2, tmp_path):
    f2.save(str(tmp_path / "d"))
    assert json.loads((tmp_path / "d" / "metadata.json").read_text())["provenance"] == f2.provenance
    assert pm.read_f2(str(tmp_path / "d")).provenance == f2.provenance
    assert pickle.loads(pickle.dumps(f2)).provenance == f2.provenance


def test_a_subset_keeps_the_provenance(f2):
    assert f2.subset(["French", "Yoruba"]).provenance == f2.provenance


def test_blocks_not_extracted_have_none():
    af = pm.read_genotypes(SMALL).allele_frequencies(["French", "Yoruba"])
    assert pm.f2_from_frequencies(af).provenance is None


def popmix(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["popmix", *args], capture_output=True, text=True)


def test_a_result_file_gets_a_provenance_sidecar(tmp_path):
    out = tmp_path / "f4.tsv"
    r = popmix("qpdstat", "-i", SMALL, "--pop1", "French", "--pop2", "Yoruba", "--pop3", "Simulation.1",
               "--pop4", "Simulation", "-o", str(out))
    assert r.returncode == 0, r.stderr
    p = json.loads((tmp_path / "f4.tsv.provenance.json").read_text())
    assert p["popmix_version"] == pm.__version__
    assert p["command"] == "qpdstat"
    assert p["argv"][:2] == ["qpdstat", "-i"]
    assert p["parameters"]["input_prefix"] == SMALL
    assert p["parameters"]["output"] == str(out)
    assert p["inputs"][f"{SMALL}.geno"]["size"] == Path(f"{SMALL}.geno").stat().st_size


def test_an_output_prefix_gets_one_sidecar(tmp_path):
    r = popmix("pca", "-i", SMALL, "--pca-pops", "French,Yoruba", "--components", "2", "-o", str(tmp_path / "p"))
    assert r.returncode == 0, r.stderr
    p = json.loads((tmp_path / "p.provenance.json").read_text())
    assert p["parameters"]["outlier_sigma"] == 6.0
    assert p["parameters"]["pca_pops"] == ["French", "Yoruba"]


def test_an_f2_directory_gets_its_sidecar_inside(tmp_path):
    r = popmix("extract-f2", "-i", SMALL, "--populations", "French,Yoruba", "-o", str(tmp_path / "f2"))
    assert r.returncode == 0, r.stderr
    assert json.loads((tmp_path / "f2" / "provenance.json").read_text())["command"] == "extract-f2"
    assert pm.read_f2(str(tmp_path / "f2")).provenance["parameters"]["populations"] == ["French", "Yoruba"]


def test_output_on_stdout_writes_no_sidecar(tmp_path):
    r = subprocess.run(["popmix", "qpdstat", "-i", SMALL, "--pop1", "French", "--pop2", "Yoruba", "--pop3",
                        "Simulation.1", "--pop4", "Simulation"], capture_output=True, text=True, cwd=tmp_path)
    assert r.returncode == 0, r.stderr
    assert list(tmp_path.iterdir()) == []


def test_a_failed_command_writes_no_sidecar(tmp_path):
    r = popmix("qpdstat", "-i", SMALL, "--pop1", "Nowhere", "--pop2", "Yoruba", "--pop3", "Simulation.1",
               "--pop4", "Simulation", "-o", str(tmp_path / "f4.tsv"))
    assert r.returncode != 0
    assert not (tmp_path / "f4.tsv.provenance.json").exists()


def test_every_command_records_provenance():
    assert [name for name, c in cli.commands.items() if not getattr(c.callback, "records_provenance", False)] == []
