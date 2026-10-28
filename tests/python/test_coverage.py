"""Coverage-aware handling of low-coverage (ancient) samples (TASKS.md H10).

``popmix samples`` / ``pm.sample_coverage`` give each sample's SNPs called,
heterozygosity and the ploidy each rule finds; PCA warns when low-coverage
samples build the axes, and KING-robust when pairs share few SNPs.
"""

from __future__ import annotations

import json
import subprocess
import warnings
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")


def caught(f) -> list[str]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        f()
    return [str(x.message) for x in w if issubclass(x.category, pm.PopmixWarning)]


@pytest.fixture(scope="module")
def ancient(tmp_path_factory) -> str:
    """The small data with French HGDP00521 missing 60% of its calls and
    Yoruba's first sample pseudo-haploid (heterozygous calls read as 2)."""
    prefix = str(tmp_path_factory.mktemp("ancient") / "a")
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", prefix, "-of", "eigenstrat"], check=True,
                   capture_output=True)
    ind = [line.split() for line in Path(f"{prefix}.ind").read_text().splitlines()]
    gappy = next(k for k, r in enumerate(ind) if r[0] == "HGDP00521")
    haploid = next(k for k, r in enumerate(ind) if r[2] == "Yoruba")
    lines = Path(f"{prefix}.geno").read_text().splitlines()
    out = []
    for n, line in enumerate(lines):
        row = list(line)
        if n % 10 < 6:
            row[gappy] = "9"
        if row[haploid] == "1":
            row[haploid] = "2"
        out.append("".join(row))
    Path(f"{prefix}.geno").write_text("\n".join(out) + "\n")
    return prefix


def test_the_table_counts_the_calls(ancient):
    df = pm.sample_coverage(ancient).set_index("id")
    geno = pm.read_genotypes(ancient)
    data = geno.to_numpy()
    for k, s in enumerate(geno.samples[:5]):
        called = int((data[:, k] < 3).sum())
        assert df.loc[s.id, "snps_called"] == called
        assert df.loc[s.id, "heterozygosity"] == pytest.approx((data[:, k] == 1).sum() / called, rel=1e-12)
    assert df.loc["HGDP00521", "called_fraction"] == pytest.approx(0.4, abs=0.01)


def test_the_table_gives_each_rules_ploidy(ancient):
    df = pm.sample_coverage(ancient, ["Yoruba"])
    first = df.iloc[0]
    assert (first["heterozygosity"], first["ploidy"], first["ploidy_admixtools2"]) == (0.0, 1, 1)
    assert set(df["ploidy"].iloc[1:]) == {2}


def test_low_coverage_pca_samples_warn(ancient):
    messages = caught(lambda: pm.smartpca(ancient, ["French", "Yoruba"], n_components=2, outlier_iterations=0))
    assert messages == ["1 PCA sample(s) miss more than half of the SNPs (up to 60%): the axes are best built "
                        "from well-covered samples, with low-coverage ones projected (--project-pops, with "
                        "--lsqproject)"]


def test_projecting_them_does_not(ancient, tmp_path):
    keep = tmp_path / "keep.ind"
    keep.write_text("".join(line + "\n" for line in Path(f"{ancient}.ind").read_text().splitlines()
                            if "HGDP00521" not in line))
    sub = str(tmp_path / "sub")
    subprocess.run(["popmix", "convertf", "-i", ancient, "-o", sub, "-of", "eigenstrat", "--keep-inds", str(keep)],
                   check=True, capture_output=True)
    assert caught(lambda: pm.smartpca(sub, ["French", "Yoruba"], n_components=2, outlier_iterations=0)) == []


def test_king_warns_on_pairs_sharing_few_snps(tmp_path):
    sub = str(tmp_path / "few")
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", sub, "-of", "packedancestrymap", "--chromosomes", "22"],
                   check=True, capture_output=True)
    r = subprocess.run(["popmix", "kinship", "-i", sub, "--populations", "French"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "share fewer than 10000 SNPs" in r.stderr


def test_king_with_enough_snps_does_not(tmp_path):
    r = subprocess.run(["popmix", "kinship", "-i", SMALL, "--populations", "French"], capture_output=True, text=True)
    assert "share fewer than" not in r.stderr


def test_the_samples_command(ancient):
    r = subprocess.run(["popmix", "samples", "-i", ancient, "--populations", "Yoruba", "--format", "json"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = json.loads(r.stdout)
    assert rows[0]["ploidy"] == 1
    assert len(rows) == 21
    assert rows[1]["called_fraction"] == pytest.approx(1.0, abs=0.01)
