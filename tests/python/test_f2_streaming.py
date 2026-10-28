"""f2 directories written a chromosome at a time and read in part (TASKS.md F4).

``extract_f2_to_dir`` (and ``popmix extract-f2``) hold one chromosome's
frequencies and blocks at a time, so the memory does not grow with every
block of every pair; the directory must be the same bytes as
``extract_f2(...).save``. ``read_f2(dir, populations=...)`` reads only
those populations' pairs and must equal reading everything and
subsetting.
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
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]
FILES = ["metadata.json", "block_sums.npy", "block_raw_sums.npy", "block_weights.npy", "block_het_sums.npy",
         "block_het_weights.npy", "block_fst_num.npy", "block_fst_den.npy"]


def same_directory(a: Path, b: Path) -> None:
    for name in FILES:
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


def both(tmp_path: Path, prefix: str, pops: list[str], **kw) -> tuple[Path, Path]:
    streamed, whole = tmp_path / "streamed", tmp_path / "whole"
    pm.extract_f2_to_dir(prefix, str(streamed), pops, **kw)
    pm.extract_f2(prefix, pops, **kw).save(str(whole), source=prefix)
    return streamed, whole


@pytest.mark.parametrize("kw", [{}, {"maxmiss": 0.5}, {"transitions": False}, {"block_size": 0.02}])
def test_the_streamed_directory_is_the_saved_one(tmp_path, kw):
    same_directory(*both(tmp_path, SMALL, POPS, **kw))


@pytest.mark.parametrize("fmt", ["eigenstrat", "ped"])
def test_other_formats_too(tmp_path, fmt):
    prefix = str(tmp_path / "data")
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", prefix, "-of", fmt], check=True, capture_output=True)
    same_directory(*both(tmp_path, prefix, POPS))


def test_a_second_pass_when_the_kept_snps_have_no_map(tmp_path):
    # every SNP at 0 cM but the first ten, then those ten left out: all SNPs say "a map", the kept ones
    # "no map", so the guessed blocks are wrong and a second pass in base pairs gives the whole's
    prefix = str(tmp_path / "data")
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", prefix, "-of", "eigenstrat"], check=True,
                   capture_output=True)
    rows = [line.split() for line in Path(f"{prefix}.snp").read_text().splitlines()]
    for r in rows[10:]:
        r[2] = "0"
    Path(f"{prefix}.snp").write_text("".join(" ".join(r) + "\n" for r in rows))
    streamed, whole = both(tmp_path, prefix, POPS, excludesnps=[r[0] for r in rows[:10]])
    same_directory(streamed, whole)
    lengths = json.loads((whole / "metadata.json").read_text())["block_lengths"]
    assert len(lengths) == len(pm.extract_f2(prefix, POPS, excludesnps=[r[0] for r in rows[:10]]).block_lengths)


def test_reading_some_populations_is_subsetting(tmp_path):
    d = tmp_path / "f2"
    pm.extract_f2_to_dir(SMALL, str(d), POPS)
    part = pm.read_f2(str(d), populations=["Simulation.1", "French"])
    full = pm.read_f2(str(d)).subset(["Simulation.1", "French"])
    assert part.populations == ["Simulation.1", "French"]
    for a, b in [("Simulation.1", "French"), ("French", "French")]:
        np.testing.assert_array_equal(part.f2_per_block(a, b), full.f2_per_block(a, b))
    assert part.fst("Simulation.1", "French") == full.fst("Simulation.1", "French")
    assert part.provenance == full.provenance


def test_an_unknown_population_is_an_error(tmp_path):
    d = tmp_path / "f2"
    pm.extract_f2_to_dir(SMALL, str(d), POPS)
    with pytest.raises(pm.PopulationNotFoundError):
        pm.read_f2(str(d), populations=["French", "Nowhere"])


def test_the_cli_writes_and_reads_in_part(tmp_path):
    d = tmp_path / "f2"
    r = subprocess.run(["popmix", "extract-f2", "-i", SMALL, "--populations", ",".join(POPS), "-o", str(d)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert pm.read_f2_populations(str(d)) == POPS
    r = subprocess.run(["popmix", "f4", "--f2-dir", str(d), "--pop1", "French", "--pop2", "Yoruba", "--pop3",
                        "Simulation", "--pop4", "Simulation.1", "--format", "json"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    want = pm.extract_f2(SMALL, POPS).f4("French", "Yoruba", "Simulation", "Simulation.1")
    assert json.loads(r.stdout)["est"] == pytest.approx(want["est"], rel=1e-12)
