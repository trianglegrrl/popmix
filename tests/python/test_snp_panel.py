"""The SNP panel is recorded, and mixing panels is said (TASKS.md H9).

f-statistics depend on how the SNPs were ascertained, so an f2 extraction
records a fingerprint of the SNP panel it read. popmix combines panels in
one place, ``merge``, which keeps only the SNPs both datasets share: when
they share less than 90% of the smaller panel, the merged data are on a
panel of their own, and it warns.
"""

from __future__ import annotations

import subprocess
import warnings
from pathlib import Path

import pytest

import popmix as pm

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")


@pytest.fixture(scope="module")
def eig(tmp_path_factory) -> str:
    """The small data as EIGENSTRAT text, to cut by hand."""
    prefix = str(tmp_path_factory.mktemp("eig") / "e")
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", prefix, "-of", "eigenstrat"], check=True,
                   capture_output=True)
    return prefix


def cut(eig: str, out: str, snps: slice, pops: set[str]) -> str:
    """EIGENSTRAT data of some SNPs (rows) and the samples of some populations."""
    ind = Path(f"{eig}.ind").read_text(encoding="utf-8").splitlines()
    keep = [k for k, line in enumerate(ind) if line.split()[2] in pops]
    Path(f"{out}.ind").write_text("".join(ind[k] + "\n" for k in keep), encoding="utf-8")
    Path(f"{out}.snp").write_text("".join(x + "\n" for x in Path(f"{eig}.snp").read_text(encoding="utf-8").splitlines()[snps]), encoding="utf-8")
    rows = Path(f"{eig}.geno").read_text(encoding="utf-8").splitlines()[snps]
    Path(f"{out}.geno").write_text("".join("".join(r[k] for k in keep) + "\n" for r in rows), encoding="utf-8")
    return out


def test_an_extraction_records_its_panel():
    panel = pm.extract_f2(SMALL, ["French", "Yoruba"]).provenance["snp_panel"]
    assert panel["snps"] == len(Path(f"{SMALL}.snp").read_text(encoding="utf-8").splitlines())
    assert len(panel["fingerprint"]) == 16


def test_the_fingerprint_is_the_panels_not_the_populations():
    a = pm.extract_f2(SMALL, ["French", "Yoruba"]).provenance["snp_panel"]
    b = pm.extract_f2(SMALL, ["Simulation.1", "French"]).provenance["snp_panel"]
    assert a == b


def test_another_panel_has_another_fingerprint(eig, tmp_path):
    n = len(Path(f"{eig}.snp").read_text(encoding="utf-8").splitlines())
    fewer = cut(eig, str(tmp_path / "fewer"), slice(0, n - 1), {"French", "Yoruba"})
    a = pm.extract_f2(SMALL, ["French", "Yoruba"]).provenance["snp_panel"]
    b = pm.extract_f2(fewer, ["French", "Yoruba"]).provenance["snp_panel"]
    assert b["snps"] == a["snps"] - 1
    assert b["fingerprint"] != a["fingerprint"]


def merge(a: str, b: str, out: str) -> subprocess.CompletedProcess:
    return subprocess.run(["popmix", "merge", "-a", a, "-b", b, "-o", out, "-of", "eigenstrat"],
                          capture_output=True, text=True)


def test_merging_partly_overlapping_panels_warns(eig, tmp_path):
    n = len(Path(f"{eig}.snp").read_text(encoding="utf-8").splitlines())
    a = cut(eig, str(tmp_path / "a"), slice(0, n * 7 // 10), {"French"})
    b = cut(eig, str(tmp_path / "b"), slice(n * 3 // 10, n), {"Yoruba"})
    shared = n * 7 // 10 - n * 3 // 10
    smaller = min(n * 7 // 10, n - n * 3 // 10)
    r = merge(a, b, str(tmp_path / "m"))
    assert r.returncode == 0, r.stderr
    assert (f"Warning: file A and file B share {shared} SNPs, {shared / smaller:.0%} of the smaller panel "
            f"({smaller}): the merged data keep only these, a panel of their own, so compare f-statistics "
            "only among populations taken from the merged data") in " ".join(r.stderr.split())


def test_merging_one_panel_does_not_warn(eig, tmp_path):
    a = cut(eig, str(tmp_path / "a"), slice(None), {"French"})
    b = cut(eig, str(tmp_path / "b"), slice(None), {"Yoruba"})
    r = merge(a, b, str(tmp_path / "m"))
    assert r.returncode == 0, r.stderr
    assert "Warning" not in r.stderr


def test_the_python_merge_warns_too(eig, tmp_path):
    n = len(Path(f"{eig}.snp").read_text(encoding="utf-8").splitlines())
    a = cut(eig, str(tmp_path / "a"), slice(0, n // 2), {"French"})
    b = cut(eig, str(tmp_path / "b"), slice(n // 4, n), {"Yoruba"})
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        pm.merge(a, b, str(tmp_path / "m"), "eigenstrat")
    assert [x.category for x in w] == [pm.PopmixWarning]
