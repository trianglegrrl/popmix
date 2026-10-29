"""A command's help says what the command prints (TASKS.md J5)."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from popmix.cli import cli


def _help(*args: str) -> str:
    result = CliRunner().invoke(cli, [*args, "--help"])
    assert result.exit_code == 0
    return result.output


def test_f4_help_says_f4_not_d():
    """f4 prints f4(A, B; C, D), not D: the help names D only to point at qpdstat."""
    text = _help("f4")
    summary = text.split("\n\n")[1]
    assert "f4" in summary
    assert "D-statistic" not in summary
    assert "ABBA" not in summary
    assert "qpdstat --d-mode" in " ".join(text.split())


def test_the_command_list_calls_f4_f4():
    line = next(line for line in _help().splitlines() if line.split()[:1] == ["f4"])
    assert "D-statistic" not in line


def test_the_readme_calls_f4_f4():
    readme = (Path(__file__).parents[2] / "README.md").read_text(encoding="utf-8")
    f4_lines = [line for line in readme.splitlines() if "`f4`" in line or line.startswith("Compute f4")]
    assert f4_lines
    for line in f4_lines:
        assert "D-statistic" not in line
