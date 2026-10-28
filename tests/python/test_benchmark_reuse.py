"""Rerunning popmix's side of the benchmarks against earlier reference runs."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[2] / "benchmarks"))

import run  # noqa: E402
import scaling  # noqa: E402


def test_a_case_reuses_the_earlier_reference_timing_and_its_outputs(tmp_path):
    ref_dir = tmp_path / "reference"
    ref_dir.mkdir()
    (ref_dir / "out.json").write_text("{}")
    earlier = {"date": "2026-10-25", "cases": {"qpadm": {"reference_runs": {"wall_s": {"median": 21.5}}}}}
    got = run.reused_reference(earlier, "qpadm", ref_dir)
    assert got == {"wall_s": {"median": 21.5}, "reused_from": "2026-10-25"}


def test_reuse_needs_the_earlier_case_and_its_outputs(tmp_path):
    earlier = {"date": "2026-10-25", "cases": {"qpadm": {"reference_runs": {}}}}
    with pytest.raises(SystemExit, match="pca"):
        run.reused_reference(earlier, "pca", tmp_path)
    with pytest.raises(SystemExit, match="no reference outputs"):
        run.reused_reference(earlier, "qpadm", tmp_path / "missing")


def test_a_scaling_point_reuses_its_reference(tmp_path):
    earlier = {"date": "2026-10-25", "curves": {"samples": {"pca:1000": {"reference": {"0": {"wall_s": 1}}}}}}
    (tmp_path / "reference-t0").mkdir()
    (tmp_path / "reference-t0" / "out.eval").write_text("1\n")
    got = scaling.reused_reference(earlier, "samples", "pca:1000", tmp_path)
    assert got == {"0": {"wall_s": 1, "reused_from": "2026-10-25"}}
    assert json.dumps(got)


def test_outputs_that_differ_only_in_provenance_digest_the_same(tmp_path):
    # the provenance records the thread count, path and time of each run
    a, b = tmp_path / "a", tmp_path / "b"
    for d, threads in ((a, 1), (b, 4)):
        d.mkdir()
        (d / "out.geno").write_bytes(b"GENO")
        (d / "out.provenance.json").write_text(json.dumps({"threads": threads}))
    assert run._digest(a) == run._digest(b)


def test_outputs_with_different_data_digest_differently(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for d, body in ((a, b"GENO1"), (b, b"GENO2")):
        d.mkdir()
        (d / "out.geno").write_bytes(body)
    assert run._digest(a) != run._digest(b)
