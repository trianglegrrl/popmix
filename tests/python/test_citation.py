"""CITATION.cff: how to cite popmix, read by GitHub and by Zenodo for each release."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
CITATION = yaml.safe_load((ROOT / "CITATION.cff").read_text())


def test_the_citation_uses_cff_1_2():
    assert CITATION["cff-version"] == "1.2.0"


def test_the_citation_version_is_the_package_version():
    pyproject = re.search(r'(?m)^version = "([^"]+)"', (ROOT / "pyproject.toml").read_text()).group(1)
    assert CITATION["version"] == pyproject


def test_the_citation_licence_is_mit_and_bsd_3():
    assert CITATION["license"] == ["MIT", "BSD-3-Clause"]


def test_the_author_is_named_with_her_orcid():
    assert CITATION["authors"] == [{"family-names": "Hardie", "given-names": "Alaina",
                                    "orcid": "https://orcid.org/0000-0002-1194-8003"}]


def _references(kind: str) -> list[dict]:
    return [ref for ref in CITATION["references"] if ref["type"] == kind]


def test_every_cited_article_has_a_doi():
    assert [ref["title"] for ref in _references("article") if not ref.get("doi")] == []


def test_every_cited_tool_says_where_to_get_it():
    assert [ref["title"] for ref in _references("software")
            if not (ref.get("url") or ref.get("repository-code"))] == []


@pytest.mark.parametrize("tool", ["admixtools2", "ADMIXTOOLS", "EIGENSOFT", "DATES", "ALDER", "ADMIXTURE",
                                  "KING", "READv2", "PLINK 2.0", "PLINK 1.9", "dadi", "BCFtools", "msprime",
                                  "R:", "L-BFGS-B", "mt19937ar"])
def test_each_reference_tool_is_cited(tool):
    assert [ref["title"] for ref in _references("software") if ref["title"].startswith(tool)] != []


def test_the_aadr_release_the_goldens_use_is_cited():
    assert [ref["doi"] for ref in _references("data")] == ["10.7910/DVN/FFIDCW"]
