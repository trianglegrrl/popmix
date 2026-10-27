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


# popmix's own exception classes (TASKS.md I2): each is a popmix.PopmixError
# and also the built-in type the error was before, so `except ValueError`
# code keeps working.
SMALL = "tests/golden/dates/family_chr20_22"
BUILTIN_OF = {
    "InputError": ValueError,
    "PopulationNotFoundError": ValueError,
    "MissingFileError": FileNotFoundError,
    "FileError": OSError,
    "ComputationError": RuntimeError,
    "InternalError": RuntimeError,
}


@pytest.mark.parametrize("name", BUILTIN_OF)
def test_each_error_class_is_a_popmix_error_and_its_builtin(name):
    cls = getattr(pm, name)
    assert issubclass(cls, pm.PopmixError)
    assert issubclass(cls, BUILTIN_OF[name])
    assert name in pm.__all__
    assert "PopmixError" in pm.__all__
    assert cls.__module__ == "popmix.errors"


def test_errors_from_the_rust_core_are_popmix_classes(tmp_path):
    with pytest.raises(pm.MissingFileError):
        pm.read_eigenstrat(str(tmp_path / "nothing"))
    prefix = tmp_path / "bad"
    Path(f"{prefix}.snp").write_text("rs1 1 0.0\n")
    Path(f"{prefix}.ind").write_text("a U P\n")
    Path(f"{prefix}.geno").write_text("0\n")
    with pytest.raises(pm.InputError, match="line 1"):
        pm.read_eigenstrat(str(prefix))
    with pytest.raises(pm.InputError, match="thread count"):
        pm.init_threads(0)


@pytest.mark.skipif(not Path(f"{SMALL}.geno").exists(), reason="test data not available")
def test_an_unknown_population_is_population_not_found():
    with pytest.raises(pm.PopulationNotFoundError, match="Nowhere"):
        pm.extract_f2(SMALL, ["French", "Nowhere"])
    f2 = pm.extract_f2(SMALL, ["French", "Yoruba"])
    with pytest.raises(pm.PopulationNotFoundError, match="Nowhere"):
        f2.f2("French", "Nowhere")


def test_the_cli_keeps_the_traceback_of_an_internal_error(monkeypatch):
    from click.testing import CliRunner

    import popmix.cli as cli_module

    def bug(*args, **kwargs):
        raise pm.InternalError("Internal error (please report it): x")

    monkeypatch.setattr(cli_module.pm, "dataset_info", bug)
    run = CliRunner().invoke(cli_module.cli, ["info", "x"])
    assert isinstance(run.exception, pm.InternalError)
