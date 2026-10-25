"""Errors keep their kind across the bindings.

A missing file is a FileNotFoundError, a malformed one or a bad argument a
ValueError; RuntimeError is left for numerical failures (a singular system).
"""

from pathlib import Path

import pytest

import popmix as pm


def test_a_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        pm.read_eigenstrat(str(tmp_path / "nothing"))
    with pytest.raises(FileNotFoundError):
        pm.read_plink(str(tmp_path / "nothing"))


def test_a_malformed_file(tmp_path):
    prefix = tmp_path / "bad"
    Path(f"{prefix}.snp").write_text("rs1 1 0.0\n")
    Path(f"{prefix}.ind").write_text("a U P\n")
    Path(f"{prefix}.geno").write_text("0\n")
    with pytest.raises(ValueError, match="line 1"):
        pm.read_eigenstrat(str(prefix))


def test_a_bad_thread_count():
    with pytest.raises(ValueError, match="thread count"):
        pm.init_threads(0)
