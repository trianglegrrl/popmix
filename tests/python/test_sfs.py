"""SFS API and the sfs command, on the 50-sample family data.

The spectra are checked against admixtools2 and dadi in tests/golden/test_sfs.py.
"""

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_projection_keeps_every_snp_with_enough_data():
    full = pm.sfs(SMALL, ["Yoruba"])
    small = pm.sfs(SMALL, ["Yoruba"], [10])
    assert full.sample_sizes == [42]
    assert small.data.sum() > full.data.sum()
    # every SNP with at least 10 Yoruba alleles called
    geno = pm.read_genotypes(SMALL)
    yoruba = [k for k, smp in enumerate(geno.samples) if smp.population == "Yoruba"]
    called = (geno.to_numpy()[:, yoruba] < 3).sum(axis=1)
    assert small.data.sum() == pytest.approx(int((2 * called >= 10).sum()))


def test_folded_spectrum_counts_the_same_snps():
    u = pm.sfs(SMALL, ["French", "Yoruba"], [8, 6])
    f = pm.sfs(SMALL, ["French", "Yoruba"], [8, 6], folded=True)
    assert f.data.sum() == pytest.approx(u.data.sum())
    assert f.mask.sum() > u.mask.sum()


def test_summary_needs_one_unfolded_population():
    with pytest.raises(ValueError, match="one-population"):
        pm.sfs(SMALL, ["French", "Yoruba"], [4, 4]).pi()
    with pytest.raises(ValueError, match="unfolded"):
        pm.sfs(SMALL, ["Yoruba"], [10], folded=True).tajima_d()


def test_joint_sfs_table_lists_every_combination():
    df = pm.joint_sfs(SMALL, ["Simulation", "Simulation.1"])
    assert len(df) == 3 * 9
    assert df[["Simulation", "Simulation.1"]].iloc[1].tolist() == [0, 1]
    assert df["n"].sum() == df.attrs["in_sfs"]


def test_unknown_population_is_an_error():
    with pytest.raises(ValueError, match="not in"):
        pm.sfs(SMALL, ["Nowhere"])


def test_sfs_command(tmp_path):
    out = tmp_path / "sfs.tsv"
    r = subprocess.run(["popmix", "sfs", "--input", SMALL, "--populations", "Yoruba", "--projection", "20",
                        "--output", str(out)], check=True, capture_output=True, text=True)
    table = pd.read_csv(out, sep="\t")
    api = pm.sfs(SMALL, ["Yoruba"], [20])
    assert table["Yoruba"].tolist() == list(range(21))
    np.testing.assert_allclose(table["n"], api.data, rtol=1e-11)
    assert f"tajima_d\t{api.tajima_d():.6g}" in r.stderr


def test_sfs_command_joint_mode(tmp_path):
    out = tmp_path / "joint.tsv"
    subprocess.run(["popmix", "sfs", "--input", SMALL, "--populations", "Simulation,Simulation.1", "--joint",
                    "--output", str(out)], check=True, capture_output=True, text=True)
    table = pd.read_csv(out, sep="\t")
    assert table["n"].tolist() == pm.joint_sfs(SMALL, ["Simulation", "Simulation.1"])["n"].tolist()
