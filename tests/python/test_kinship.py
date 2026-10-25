"""Kinship API and the kinship command, on the 54-sample DATES data.

The estimates are checked against PLINK 2 (KING-robust) and READv2 in
tests/golden/test_king_plink2.py and tests/golden/test_read_readv2.py.
"""

import os
import subprocess

import pandas as pd
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_king_on_all_populations():
    df = pm.king(SMALL)
    assert len(df) == 54 * 53 // 2
    assert list(df.columns) == ["id1", "id2", "nsnp", "hethet", "ibs0", "het1_hom2", "het2_hom1", "kinship"]


def test_read_normalizes_by_the_median_pair():
    df = pm.read_kinship(SMALL, ["French", "Yoruba"])
    assert df["normalized_p0"].median() == pytest.approx(1.0)
    assert df.attrs["norm_value"] == pytest.approx(df["p0"].median())


def test_read_with_a_given_normalization():
    df = pm.read_kinship(SMALL, ["French"], norm_value=0.25)
    assert df.attrs["norm_value"] == 0.25
    assert (df["normalized_p0"] == df["p0"] / 0.25).all()


def test_read_rejects_unknown_options():
    with pytest.raises(ValueError, match="norm must be"):
        pm.read_kinship(SMALL, ["French"], norm="mode")
    with pytest.raises(ValueError, match="unknown compat"):
        pm.read_kinship(SMALL, ["French"], compat="readv1")


@pytest.mark.parametrize("method", ["king", "read"])
def test_kinship_command_matches_the_api(tmp_path, method):
    out = tmp_path / "kin.tsv"
    subprocess.run(["popmix", "kinship", "--input", SMALL, "--populations", "French,Yoruba", "--method", method,
                    "--output", str(out)], check=True, capture_output=True, text=True)
    got = pd.read_csv(out, sep="\t")
    api = (pm.king if method == "king" else pm.read_kinship)(SMALL, ["French", "Yoruba"])
    assert got["id1"].tolist() == api["id1"].tolist()
    assert got["kinship"].to_numpy() == pytest.approx(api["kinship"].to_numpy(), rel=1e-11)
