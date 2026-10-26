"""Start-up cost: pandas loads only when something needs it (TASKS.md F15)."""

import subprocess
import sys

import popmix as pm


def _modules_after(code: str) -> set[str]:
    out = subprocess.run([sys.executable, "-c", code + "\nimport sys\nprint(' '.join(sys.modules))"],
                         check=True, capture_output=True, text=True).stdout
    return set(out.split())


def test_importing_popmix_does_not_load_pandas():
    assert "pandas" not in _modules_after("import popmix")


def test_cli_help_does_not_load_pandas():
    code = ("from click.testing import CliRunner\nfrom popmix.cli import cli\n"
            "assert CliRunner().invoke(cli, ['--help']).exit_code == 0")
    assert "pandas" not in _modules_after(code)


def test_every_public_name_still_resolves():
    for name in pm.__all__:
        assert getattr(pm, name) is not None, name
    from popmix import qpadm_multi, smartpca  # noqa: F401


def test_functions_named_like_their_module_stay_functions():
    import popmix.admixture
    import popmix.rolloff
    import popmix.sfs  # noqa: F401

    assert callable(pm.admixture) and callable(pm.rolloff) and callable(pm.sfs)
    assert pm.admixture.__name__ == "admixture"
