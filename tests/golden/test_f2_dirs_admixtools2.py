"""Golden master tests: admixtools2 f2 directories.

tests/golden/admixtools2/v2.0.10/dg_6pop/f2_dir is admixtools2
extract_f2() of the dg_6pop populations (generate_f2_dir.R). popmix must
read it to admixtools2's f2 blocks (f2.json, exactly), which match popmix's
own (rtol 1e-10, counts exact); and admixtools2 must read the directory
popmix writes back to the same values (with --run-slow, if R and admixtools2
are installed).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

import popmix as pm

from .reference_paths import local

REF_DIR = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop"
F2_DIR = REF_DIR / "f2_dir"
RSCRIPT = "/usr/bin/Rscript"

pytestmark = pytest.mark.skipif(not F2_DIR.exists(), reason="f2_dir not generated")


def _ref() -> dict:
    return json.loads((REF_DIR / "f2.json").read_text())


@pytest.fixture(scope="module")
def theirs():
    return pm.read_admixtools2_f2(F2_DIR)


@pytest.fixture(scope="module")
def ours():
    ref = _ref()
    if not os.path.exists(f"{local(ref['prefix'])}.geno"):
        pytest.skip(f"genotype data not available: {local(ref['prefix'])}")
    p = ref["params"]
    return pm.extract_f2(local(ref["prefix"]), ref["pops"], auto_only=p["auto_only"], block_size=p["blgsize"],
                         maxmiss=p["maxmiss"], poly_only=p["poly_only"] == "f2")


def test_reads_admixtools2_blocks_exactly(theirs):
    ref = _ref()
    assert sorted(theirs.populations) == sorted(ref["pops"])
    assert list(theirs.block_lengths) == ref["block_lengths"]
    for pair, blocks in ref["f2_blocks"].items():
        a, b = pair.split("|")
        np.testing.assert_array_equal(theirs.f2_per_block(a, b), np.array(blocks, dtype=float))


def test_matches_popmix_blocks(theirs, ours):
    pops = list(ours.populations)
    for i, a in enumerate(pops):
        for b in pops[i:]:
            np.testing.assert_allclose(theirs.f2_per_block(a, b), ours.f2_per_block(a, b), rtol=1e-10, atol=1e-15)
            # admixtools2 stores a population with itself as fully covered
            if a != b:
                assert theirs.counts_per_block(a, b) == ours.counts_per_block(a, b), (a, b)
    t, o = pm.qpadm(theirs, ["Han.DG", "Yoruba.DG"], ["Mbuti.DG", "Papuan.DG", "Chimp.REF"], "French.DG"), \
        pm.qpadm(ours, ["Han.DG", "Yoruba.DG"], ["Mbuti.DG", "Papuan.DG", "Chimp.REF"], "French.DG")
    np.testing.assert_allclose(t.weights, o.weights, rtol=1e-9)


def test_subset_and_population_selection(theirs):
    sub = pm.read_admixtools2_f2(F2_DIR, ["Yoruba.DG", "Han.DG", "Chimp.REF"])
    assert list(sub.populations) == ["Yoruba.DG", "Han.DG", "Chimp.REF"]
    np.testing.assert_array_equal(sub.f2_per_block("Han.DG", "Yoruba.DG"), theirs.f2_per_block("Han.DG", "Yoruba.DG"))
    picked = pm.extract_f2_subset(theirs, ["Han.DG", "Yoruba.DG"])
    assert list(picked.populations) == ["Han.DG", "Yoruba.DG"]
    assert picked.f2("Han.DG", "Yoruba.DG") == theirs.f2("Han.DG", "Yoruba.DG")


def test_write_then_read_round_trips(theirs, tmp_path):
    out = tmp_path / "f2"
    pm.write_admixtools2_f2(theirs, out)
    back = pm.read_admixtools2_f2(out)
    for a in theirs.populations:
        for b in theirs.populations:
            np.testing.assert_array_equal(back.f2_per_block(a, b), theirs.f2_per_block(a, b))
    with pytest.raises(FileExistsError):
        pm.write_admixtools2_f2(theirs, out)


@pytest.mark.slow
@pytest.mark.skipif(not os.path.exists(RSCRIPT), reason="R not installed")
def test_admixtools2_reads_what_popmix_writes(ours, tmp_path):
    out = tmp_path / "popmix_f2"
    pm.write_admixtools2_f2(ours, out)
    script = f"""
    suppressPackageStartupMessages(library(admixtools))
    f2 <- f2_from_precomp("{out}", verbose = FALSE)
    cat(jsonlite::toJSON(list(pops = dimnames(f2)[[1]], lengths = as.numeric(gsub("l", "", dimnames(f2)[[3]])),
                              han_yoruba = sprintf("%.17g", f2["Han.DG", "Yoruba.DG", ]))))
    """
    if shutil.which(RSCRIPT) is None:
        pytest.skip("R not installed")
    res = subprocess.run([RSCRIPT, "-e", script], capture_output=True, text=True)
    if "there is no package called" in res.stderr:
        pytest.skip("admixtools2 not installed")
    assert res.returncode == 0, res.stderr
    got = json.loads(res.stdout)
    assert sorted(got["pops"]) == sorted(ours.populations)
    assert got["lengths"] == list(ours.block_lengths)
    np.testing.assert_array_equal([float(x) for x in got["han_yoruba"]], ours.f2_per_block("Han.DG", "Yoruba.DG"))
