"""The CLI reports bad input briefly and lets bugs show (B12).

Every command used to catch every exception and print one line, so a bug
(a TypeError, say) looked like bad input and lost its traceback.
"""

import pytest
from click.testing import CliRunner

import popmix.cli as cli_module
import popmix.cli._population_checks as checks
from popmix.cli import cli


def test_bad_input_is_one_line(tmp_path):
    run = CliRunner().invoke(cli, ["info", str(tmp_path / "nothing")])
    assert run.exit_code == 1
    assert run.output.startswith("Error: ")
    assert "Traceback" not in run.output


def test_debug_shows_the_traceback(tmp_path):
    run = CliRunner().invoke(cli, ["--debug", "info", str(tmp_path / "nothing")])
    assert run.exit_code == 1
    assert isinstance(run.exception, FileNotFoundError)


def test_a_bug_is_not_reported_as_bad_input(monkeypatch, tmp_path):
    def broken(prefix):
        raise TypeError("a bug")

    monkeypatch.setattr(cli_module.pm, "dataset_info", broken)
    run = CliRunner().invoke(cli, ["info", str(tmp_path / "x")])
    assert isinstance(run.exception, TypeError)
    assert "Error: a bug" not in run.output


@pytest.mark.parametrize("handler", ["suggest", "validate"])
def test_helpers_only_swallow_expected_errors(monkeypatch, handler):
    def broken(*args, **kwargs):
        raise TypeError("a bug")

    if handler == "suggest":
        monkeypatch.setattr(cli_module.pm, "suggest_populations", broken)
        with pytest.raises(TypeError):
            checks.format_population_error("X", "prefix")
    else:
        monkeypatch.setattr(cli_module.pm, "list_populations", broken)
        with pytest.raises(TypeError):
            checks.validate_populations("prefix", ["X"])
