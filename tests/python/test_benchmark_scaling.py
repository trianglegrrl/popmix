"""The scaling benchmark's input series (benchmarks/scaling.py, TASKS.md E5)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[2] / "benchmarks"))

import scaling  # noqa: E402

pytestmark = pytest.mark.skipif(not Path(f"{scaling.AADR}.ind").exists() or not scaling.ANNO.exists(),
                                reason="AADR v66.p1 not available")


def test_population_series_has_enough_samples_per_group():
    groups = scaling.pop_groups()
    assert len(groups) == max(scaling.POP_COUNTS)
    counts: dict[str, int] = {}
    for _, _, g in scaling._ind_rows():
        counts[g] = counts.get(g, 0) + 1
    assert all(counts[g] >= scaling.PER_POP for g in groups)
    assert groups == scaling.pop_groups()  # a fixed order, so the series is nested


def test_sample_series_is_every_present_day_sample():
    ids = scaling.present_day_samples()
    assert len(ids) == len(set(ids)) == max(scaling.SAMPLE_COUNTS)
    ind = {r[0] for r in scaling._ind_rows()}
    assert set(ids) <= ind


def test_point_keys_are_unique_within_each_curve():
    seen = set()
    for p in scaling.points():
        assert (p.curve, p.key) not in seen
        seen.add((p.curve, p.key))
    assert {p.curve for p in scaling.points()} == {"populations", "snps", "samples", "threads", "full_aadr"}


def test_sample_series_is_one_population(tmp_path, monkeypatch):
    """smartpca stops at more than 100 populations, so the sample series
    labels every kept sample one population."""
    written = {}

    def fake_run(cmd, **kwargs):
        written["ind"] = (tmp_path / "samples3_1pop.in.ind").read_text()

    monkeypatch.setattr(scaling, "SCALING", tmp_path)
    monkeypatch.setattr(scaling.subprocess, "run", fake_run)
    scaling._samples_input(3)()
    labels = [line.split()[2] for line in written["ind"].splitlines()]
    assert labels.count(scaling.PRESENT_DAY_LABEL) == 3
    assert set(labels) == {scaling.PRESENT_DAY_LABEL, "Ignore"}
