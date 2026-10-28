"""Choosing and relabelling samples and SNPs (TASKS.md H2, P3).

convertf keeps or leaves out samples by ID, gives them new populations,
chooses them by population, and keeps SNPs by ID, chromosome and region;
each output is checked against the input cut by hand, on the streamed path
(SNP-major output) and the in-memory one (PED, TGENO). The f2 commands take
the same sample options.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")
SOURCE = pm.read_genotypes(SMALL)
DATA = SOURCE.to_numpy()
IDS = [s.id for s in SOURCE.samples]
POPS = [s.population for s in SOURCE.samples]
SNP_IDS = [s.id for s in SOURCE.snps]


def popmix(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["popmix", *args], capture_output=True, text=True)


def convert(tmp_path: Path, *args: str, fmt: str = "packedancestrymap") -> pm.GenotypeMatrix:
    out = str(tmp_path / f"out_{fmt}")
    r = popmix("convertf", "-i", SMALL, "-o", out, "-of", fmt, *args)
    assert r.returncode == 0, r.stderr
    return pm.read_genotypes(out)


def file_of(tmp_path: Path, name: str, lines: list[str]) -> str:
    path = tmp_path / name
    path.write_text("".join(x + "\n" for x in lines))
    return str(path)


@pytest.mark.parametrize("fmt", ["packedancestrymap", "eigenstrat", "ped", "tgeno"])
def test_populations_and_a_region(tmp_path, fmt):
    g = convert(tmp_path, "--populations", "Yoruba,French", "--region", "21:15000000-30000000", fmt=fmt)
    cols = [k for k, p in enumerate(POPS) if p in ("French", "Yoruba")]
    rows = [k for k, s in enumerate(SOURCE.snps) if s.chromosome == "21" and 15_000_000 <= s.position_bp <= 30_000_000]
    assert [s.id for s in g.samples] == [IDS[k] for k in cols]
    assert [s.id for s in g.snps] == [SNP_IDS[k] for k in rows]
    np.testing.assert_array_equal(g.to_numpy(), DATA[np.ix_(rows, cols)])


def test_samples_by_id_then_relabelled(tmp_path):
    keep = file_of(tmp_path, "keep.txt", IDS[:30])
    exclude = file_of(tmp_path, "exclude.txt", IDS[:2])
    relabel = file_of(tmp_path, "relabel.txt", [f"{IDS[2]} Mine"])
    rename = file_of(tmp_path, "rename.txt", ["Mine Ours", "French Fr"])
    g = convert(tmp_path, "--keep-inds", keep, "--exclude-inds", exclude, "--relabel", relabel,
                "--relabel-pops", rename)
    assert [s.id for s in g.samples] == IDS[2:30]
    want = ["Ours"] + ["Fr" if p == "French" else p for p in POPS[3:30]]
    assert [s.population for s in g.samples] == want
    np.testing.assert_array_equal(g.to_numpy(), DATA[:, 2:30])


def test_snps_by_id_and_chromosome(tmp_path):
    keep = file_of(tmp_path, "keep.snp", SNP_IDS[::3])
    exclude = file_of(tmp_path, "exclude.txt", SNP_IDS[:30])
    g = convert(tmp_path, "--keep-snps", keep, "--exclude-snps", exclude, "--chromosomes", "20,22")
    rows = [k for k, s in enumerate(SOURCE.snps) if k % 3 == 0 and k >= 30 and s.chromosome in ("20", "22")]
    assert [s.id for s in g.snps] == [SNP_IDS[k] for k in rows]
    np.testing.assert_array_equal(g.to_numpy(), DATA[rows])


def test_an_unknown_sample_is_an_error(tmp_path):
    keep = file_of(tmp_path, "keep.txt", ["nobody", IDS[0]])
    r = popmix("convertf", "-i", SMALL, "-o", str(tmp_path / "o"), "-of", "eigenstrat", "--keep-inds", keep)
    assert r.returncode == 1
    assert "1 not in the dataset: nobody" in r.stderr


def test_a_reversed_region_is_a_usage_error(tmp_path):
    r = popmix("convertf", "-i", SMALL, "-o", str(tmp_path / "o"), "-of", "eigenstrat", "--region", "20:9-1")
    assert r.returncode == 2


def test_select_samples_gives_extract_f2s_inds_and_pops():
    inds, pops = pm._core.select_samples(SMALL, exclude_inds=[IDS[0]], relabel_pops={"French": "Fr"})
    assert inds == IDS[1:]
    assert pops == ["Fr" if p == "French" else p for p in POPS[1:]]


def test_f2_commands_take_a_relabelling(tmp_path):
    french = [i for i, p in zip(IDS, POPS) if p == "French"]
    relabel = file_of(tmp_path, "relabel.txt", [f"{i} FrenchB" for i in french[::2]])
    r = popmix("f2", "--pop1", "French", "--pop2", "FrenchB", "-i", SMALL, "--relabel", relabel, "--format", "json")
    assert r.returncode == 0, r.stderr
    labels = dict.fromkeys(french[::2], "FrenchB")
    want = pm.extract_f2(SMALL, ["French", "FrenchB"], inds=IDS, pops=[labels.get(i, p) for i, p in zip(IDS, POPS)])
    assert json.loads(r.stdout)["est"] == pytest.approx(want.f2("French", "FrenchB")["est"], rel=1e-12)


def test_extract_f2_takes_the_selection(tmp_path):
    rename = file_of(tmp_path, "rename.txt", ["French Fr"])
    exclude = file_of(tmp_path, "exclude.ind", [line for line in Path(f"{SMALL}.ind").read_text().splitlines()
                                                if line.split()[2] == "Simulation"])
    r = popmix("extract-f2", "-i", SMALL, "--relabel-pops", rename, "--exclude-inds", exclude,
               "-o", str(tmp_path / "f2"))
    assert r.returncode == 0, r.stderr
    assert pm.read_f2(str(tmp_path / "f2")).populations == ["Fr", "Yoruba", "Simulation.1"]


def test_a_population_the_selection_leaves_out_is_an_error(tmp_path):
    rename = file_of(tmp_path, "rename.txt", ["French Fr"])
    r = popmix("f2", "--pop1", "French", "--pop2", "Yoruba", "-i", SMALL, "--relabel-pops", rename)
    assert r.returncode == 1
    assert "Not among the selected samples' populations: French" in r.stderr
