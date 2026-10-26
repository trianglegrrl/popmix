"""Agreement checks of the benchmark cases (benchmarks/cases.py)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2] / "benchmarks"))

import cases  # noqa: E402


def test_f2_mismatch_is_per_entry():
    """An exact zero in the reference is not an infinite relative error."""
    want = np.array([0.0, 1.0, 2.0])
    got = np.array([1e-16, 1.0 + 5e-11, 2.0])
    assert cases._f2_mismatches(got, want) == 0


def test_f2_mismatch_counts_entries_off_both_tolerances():
    want = np.array([0.0, 1.0, 0.94664])
    got = np.array([1e-14, 1.0 + 1e-9, 0.9])
    assert cases._f2_mismatches(got, want) == 3


def test_f2_mismatch_needs_matching_nans():
    assert cases._f2_mismatches(np.array([np.nan, 1.0]), np.array([np.nan, 1.0])) == 0
    assert cases._f2_mismatches(np.array([np.nan, 1.0]), np.array([0.5, 1.0])) == 1


def test_f2_case_reads_ploidy_as_admixtools2():
    """Diploid genomes admixtools2 takes as pseudo-haploid (DEVIATIONS.md ADMIXTOOLS2-0007) would differ."""
    cmd = cases.F2_CASE.popmix(Path("in"), Path("out"), 4)
    assert cmd[cmd.index("--pseudohaploid") + 1] == "admixtools2"


def test_seed_runs_run_each_seed_once_and_cache(tmp_path):
    calls = []

    def cmd(d: Path, seed: int) -> list[str]:
        calls.append(seed)
        return ["sh", "-c", f"echo {seed * 10} > {d / 'out.txt'}"]

    def read(d: Path) -> float:
        return float((d / "out.txt").read_text())

    assert cases._seed_runs(tmp_path, "tool", cmd, read) == {s: s * 10.0 for s in cases.SPREAD_SEEDS}
    assert cases._seed_runs(tmp_path, "tool", cmd, read) == {s: s * 10.0 for s in cases.SPREAD_SEEDS}
    assert calls == list(cases.SPREAD_SEEDS)


def test_stochastic_reference_scripts_take_the_seed():
    """qpGraph's random starts and find_graphs' search depend on the seed (E3)."""
    for script in (cases._R_QPGRAPH, cases._R_FIND_GRAPHS):
        assert "seed" in script and "a[1]" in script
    qp = cases._R_QPGRAPH
    assert "seed = as.integer(a[4])" in qp
    assert "set.seed(as.integer(a[5]))" in cases._R_FIND_GRAPHS
