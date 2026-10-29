"""The GitHub Actions release workflow builds what PyPI users need and
publishes only from a version tag."""

from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/release.yml"
# one stable-ABI (abi3) wheel per platform serves every Python from the
# oldest supported on; each is installed and run on the oldest and newest
OLDEST, NEWEST = "3.9", "3.14"


@pytest.fixture(scope="module")
def workflow() -> dict:
    return yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))


def _triggers(workflow: dict) -> dict:
    # YAML 1.1 reads the bare key `on` as True
    return workflow.get("on", workflow.get(True))


def _runners(job: dict) -> list[str]:
    include = job.get("strategy", {}).get("matrix", {}).get("include")
    return [m["runner"] for m in include] if include else [job["runs-on"]]


def _uses(job: dict, action: str) -> list[dict]:
    return [s for s in job["steps"] if s.get("uses", "").startswith(action)]


def test_runs_on_version_tags_and_by_hand(workflow: dict) -> None:
    on = _triggers(workflow)
    assert on["push"]["tags"] == ["v*"]
    assert "workflow_dispatch" in on


def test_builds_one_wheel_for_each_platform(workflow: dict) -> None:
    jobs = workflow["jobs"]
    for name, runner in [("linux", "ubuntu"), ("macos", "macos"), ("windows", "windows")]:
        job = jobs[name]
        assert all(r.startswith(runner) for r in _runners(job)), name
        builds = _uses(job, "PyO3/maturin-action")
        assert builds, f"{name} does not build with maturin"
        args = builds[0]["with"]["args"]
        assert "--release" in args
        assert "-i " not in f"{args} ", f"{name}: an abi3 wheel is built once, not per interpreter"
    assert "manylinux" in _uses(jobs["linux"], "PyO3/maturin-action")[0]["with"]
    # Intel and ARM wheels: one each on Linux, one universal2 wheel on macOS
    targets = {m["target"] for m in jobs["linux"]["strategy"]["matrix"]["include"]}
    assert targets == {"x86_64", "aarch64"}
    assert _uses(jobs["macos"], "PyO3/maturin-action")[0]["with"]["target"] == (
        "universal2-apple-darwin"
    )
    assert _uses(jobs["sdist"], "PyO3/maturin-action")[0]["with"]["command"] == "sdist"


def test_the_wheel_uses_the_stable_abi_from_the_oldest_python() -> None:
    maturin = (ROOT / "pyproject.toml").read_text(encoding="utf-8").split("[tool.maturin]")[1].split("\n[")[0]
    features = next(line for line in maturin.splitlines() if line.startswith("features"))
    assert f'"pyo3/abi3-py{OLDEST.replace(".", "")}"' in features
    project = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert f'requires-python = ">={OLDEST}"' in project


CHECK = "bash .github/scripts/check-wheel.sh"
SET_VERSION = "python .github/scripts/set-version.py ${{ github.ref_name }}"
# v1.0.0 goes to PyPI; a test tag v1.0.0rc1 to TestPyPI only (TASKS.md K8)
FINAL_TAG = "startsWith(github.ref, 'refs/tags/v') && !contains(github.ref_name, 'rc')"
TEST_TAG = "startsWith(github.ref, 'refs/tags/v') && contains(github.ref_name, 'rc')"


def test_wheels_are_installed_and_run_on_the_oldest_and_newest_python(workflow: dict) -> None:
    for name in ("linux", "macos", "windows"):
        steps = workflow["jobs"][name]["steps"]
        build = next(i for i, s in enumerate(steps) if s.get("uses", "").startswith("PyO3/maturin-action"))
        checked = [steps[i - 1]["with"]["python-version"] for i, s in enumerate(steps)
                   if i > build and s.get("run", "").strip() == CHECK]
        assert checked == [OLDEST, NEWEST], name
        assert all(steps[i - 1].get("uses", "").startswith("actions/setup-python")
                   for i, s in enumerate(steps) if s.get("run", "").strip() == CHECK), name


def test_the_check_installs_the_wheel_runs_it_and_finds_the_skill() -> None:
    # the skill is package data (python/popmix/skill/): the built wheel must hold it
    script = "".join((ROOT / ".github/scripts" / name).read_text(encoding="utf-8")
                     for name in ("check-wheel.sh", "check-installed.sh"))
    for line in ("pip install popmix --no-index --find-links dist --no-deps",
                 "pip install popmix --find-links dist",
                 "popmix --version",
                 "popmix install-skill --dest",
                 "commands/qpadm.md",
                 "import popmix, popmix._core"):
        assert line in script, line


def test_publishes_from_tags_with_trusted_publishing(workflow: dict) -> None:
    publish = workflow["jobs"]["publish"]
    assert set(publish["needs"]) == {"linux", "macos", "windows", "sdist"}
    assert publish["if"] == FINAL_TAG
    assert publish["permissions"]["id-token"] == "write"
    assert publish["environment"]["name"] == "pypi"
    assert _uses(publish, "pypa/gh-action-pypi-publish")
    # no long-lived token: trusted publishing only
    assert "password" not in str(publish)


def test_every_build_sets_the_version_from_the_ref_first(workflow: dict) -> None:
    for name in ("linux", "macos", "windows", "sdist"):
        steps = workflow["jobs"][name]["steps"]
        build = next(i for i, s in enumerate(steps) if s.get("uses", "").startswith("PyO3/maturin-action"))
        runs = [s.get("run", "").strip() for s in steps[:build]]
        assert SET_VERSION in runs, name


def test_test_tags_publish_to_testpypi_only(workflow: dict) -> None:
    job = workflow["jobs"]["testpypi"]
    assert job["if"] == TEST_TAG
    assert set(job["needs"]) == {"linux", "macos", "windows", "sdist"}
    assert job["environment"]["name"] == "testpypi"
    assert job["permissions"]["id-token"] == "write"
    assert _uses(job, "pypa/gh-action-pypi-publish")[0]["with"]["repository-url"] == (
        "https://test.pypi.org/legacy/")


def test_a_test_tag_is_installed_back_from_testpypi_on_every_platform(workflow: dict) -> None:
    job = workflow["jobs"]["testpypi-install"]
    assert job["needs"] == "testpypi"
    runners = set(_runners(job))
    assert runners == {"ubuntu-latest", "ubuntu-24.04-arm", "macos-15", "windows-latest"}
    script = "\n".join(s.get("run", "") for s in job["steps"])
    assert "--index-url https://test.pypi.org/simple/" in script
    assert "bash .github/scripts/check-installed.sh" in script
