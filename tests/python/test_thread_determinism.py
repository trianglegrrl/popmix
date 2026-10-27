"""Every parallel computation gives the same bytes at any thread count (TASKS.md F8).

Sums run in fixed chunks added in a fixed order, so the thread count only
changes who adds what, never the order. The Rust tests check the f2 blocks,
the frequency pass, PCA, DATES' bins and the genotype f-statistics on their
own; this runs each analysis end to end through the bindings.
"""

from __future__ import annotations

import dataclasses
import os

import numpy as np
import pandas as pd
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
THREADS = (1, 4, 7)

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def _canon(x):
    """x as nested lists of exact values: arrays and floats by their bytes."""
    if isinstance(x, pd.DataFrame):
        return [list(x.columns), [_canon(x[c].to_numpy()) for c in x.columns], dict(x.attrs)]
    if isinstance(x, np.ndarray):
        return (x.dtype.str, x.shape, x.tobytes()) if x.dtype != object else [_canon(v) for v in x]
    if isinstance(x, float):
        return np.float64(x).tobytes()
    if isinstance(x, dict):
        return sorted((k, _canon(v)) for k, v in x.items())
    if isinstance(x, (list, tuple)):
        return [_canon(v) for v in x]
    if dataclasses.is_dataclass(x):
        return _canon(dataclasses.asdict(x))
    if hasattr(x, "__dict__") and not isinstance(x, type):
        return _canon(vars(x))
    if isinstance(x, (int, str, bytes, type(None))):
        return x
    # a Rust result object: its public fields
    return [(a, _canon(getattr(x, a))) for a in dir(x) if not a.startswith("_") and not callable(getattr(x, a))]


def _at_each_thread_count(run):
    results = []
    try:
        for t in THREADS:
            pm.init_threads(t)
            results.append(_canon(run()))
    finally:
        pm.init_threads(4)
    for t, got in zip(THREADS[1:], results[1:]):
        assert got == results[0], f"{t} threads differ from 1"


def _f2_blocks():
    f2 = pm.extract_f2(SMALL, ["French", "Yoruba", "Simulation.1"])
    return [f2.f2_per_block(a, b) for a, b in [("French", "Yoruba"), ("French", "Simulation.1")]]


ANALYSES = {
    "pca": lambda: pm.smartpca(SMALL, ["French", "Yoruba"], ["Simulation.1"]),
    "king": lambda: pm.king(SMALL, ["French", "Yoruba"]),
    "read": lambda: pm.read_kinship(SMALL, ["French"]),
    "admixture": lambda: pm.admixture(SMALL, 2, ["French", "Yoruba"]),
    "rolloff": lambda: pm.rolloff(SMALL, "Simulation.1", ["French", "Yoruba"]),
    "dates": lambda: pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"]),

    "f2_blocks": _f2_blocks,
}


@pytest.mark.parametrize("name", ANALYSES)
def test_results_do_not_depend_on_the_thread_count(name):
    _at_each_thread_count(ANALYSES[name])


@pytest.mark.parametrize("fmt", ["eigenstrat", "plink", "ped", "tgeno"])
def test_converted_files_do_not_depend_on_the_thread_count(tmp_path, fmt):
    def run():
        out = tmp_path / f"out{pm.get_threads()}"
        pm.convertf(SMALL, str(out), fmt)
        return sorted((p.suffix, p.read_bytes()) for p in tmp_path.glob(f"{out.name}.*"))
    _at_each_thread_count(run)


def test_merged_files_do_not_depend_on_the_thread_count(tmp_path):
    # the family merged with a copy of itself under other sample IDs, both
    # packed, read in place
    pm.convertf(SMALL, str(tmp_path / "a"), "packedancestrymap")
    for ext in ("geno", "snp"):
        (tmp_path / f"b.{ext}").write_bytes((tmp_path / f"a.{ext}").read_bytes())
    ind = (tmp_path / "a.ind").read_text().splitlines()
    (tmp_path / "b.ind").write_text("".join(f"copy_{line.lstrip()}\n" for line in ind))

    def run():
        out = tmp_path / f"merged{pm.get_threads()}"
        stats = pm.merge(str(tmp_path / "a"), str(tmp_path / "b"), str(out), "tgeno")
        assert stats["samples_merged"] == 2 * len(ind)
        return sorted((p.suffix, p.read_bytes()) for p in tmp_path.glob(f"{out.name}.*"))
    _at_each_thread_count(run)
