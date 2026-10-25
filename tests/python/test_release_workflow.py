"""The GitHub Actions release workflow builds what PyPI users need and
publishes only from a version tag."""

from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

WORKFLOW = Path(__file__).resolve().parents[2] / ".github/workflows/release.yml"
PYTHONS = ["3.9", "3.10", "3.11", "3.12", "3.13"]


@pytest.fixture(scope="module")
def workflow() -> dict:
    return yaml.safe_load(WORKFLOW.read_text())


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


def test_builds_wheels_for_each_platform_and_python(workflow: dict) -> None:
    jobs = workflow["jobs"]
    for name, runner in [("linux", "ubuntu"), ("macos", "macos"), ("windows", "windows")]:
        job = jobs[name]
        assert all(r.startswith(runner) for r in _runners(job)), name
        builds = _uses(job, "PyO3/maturin-action")
        assert builds, f"{name} does not build with maturin"
        args = builds[0]["with"]["args"]
        assert "--release" in args
        for python in PYTHONS:
            assert python in args, f"{name} misses Python {python}"
    assert "manylinux" in _uses(jobs["linux"], "PyO3/maturin-action")[0]["with"]
    # Intel and ARM wheels: one each on Linux, one universal2 wheel on macOS
    targets = {m["target"] for m in jobs["linux"]["strategy"]["matrix"]["include"]}
    assert targets == {"x86_64", "aarch64"}
    assert _uses(jobs["macos"], "PyO3/maturin-action")[0]["with"]["target"] == (
        "universal2-apple-darwin"
    )
    assert _uses(jobs["sdist"], "PyO3/maturin-action")[0]["with"]["command"] == "sdist"


def test_wheels_are_installed_and_run_before_publishing(workflow: dict) -> None:
    for name in ("linux", "macos", "windows"):
        runs = " ".join(s.get("run", "") for s in workflow["jobs"][name]["steps"])
        assert "pip install" in runs and "popmix --version" in runs, name


def test_publishes_from_tags_with_trusted_publishing(workflow: dict) -> None:
    publish = workflow["jobs"]["publish"]
    assert set(publish["needs"]) == {"linux", "macos", "windows", "sdist"}
    assert publish["if"] == "startsWith(github.ref, 'refs/tags/v')"
    assert publish["permissions"]["id-token"] == "write"
    assert publish["environment"]["name"] == "pypi"
    assert _uses(publish, "pypa/gh-action-pypi-publish")
    # no long-lived token: trusted publishing only
    assert "password" not in str(publish)
