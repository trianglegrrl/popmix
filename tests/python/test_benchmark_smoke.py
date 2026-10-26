"""The CI subset of the benchmarks (benchmarks/smoke.py, TASKS.md E6)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "benchmarks"))

import smoke  # noqa: E402


def test_every_command_runs_on_the_shipped_data():
    times = smoke.run(repeats=1)
    assert list(times) == ["extract-f2", "qpdstat", "pca", "kinship", "dates"]
    assert all(len(ts) == 1 and ts[0] > 0 for ts in times.values())


def test_render_is_one_row_per_command():
    table = smoke.render({"pca": [1.0, 2.0, 3.0], "dates": [0.5]})
    assert "| pca | 2.00 s | 1.00 s | 3.00 s |" in table
    assert "| dates | 0.50 s | 0.50 s | 0.50 s |" in table
