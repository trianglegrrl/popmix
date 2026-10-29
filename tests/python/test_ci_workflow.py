"""The CI workflow tests every platform a wheel is built for, and the goldens
that ship their data cannot skip there (TASKS.md K1)."""

from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/ci.yml"
RELEASE = ROOT / ".github/workflows/release.yml"


@pytest.fixture(scope="module")
def workflow() -> dict:
    return yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))


def _runners(job: dict) -> list[str]:
    include = job.get("strategy", {}).get("matrix", {}).get("include")
    return [m["runner"] for m in include] if include else [job["runs-on"]]


def _step(job: dict, name: str) -> dict:
    (step,) = [s for s in job["steps"] if s.get("name") == name]
    return step


def test_the_tests_run_on_every_platform_a_wheel_is_built_on(workflow: dict) -> None:
    release = yaml.safe_load(RELEASE.read_text(encoding="utf-8"))["jobs"]
    built = {r for job in ("linux", "macos", "windows") for r in _runners(release[job])}
    assert built <= set(_runners(workflow["jobs"]["test"]))


def test_the_universal2_macos_wheel_is_tested_on_intel_too(workflow: dict) -> None:
    assert "macos-15-intel" in _runners(workflow["jobs"]["test"])


@pytest.mark.parametrize("name", ["Rust tests", "Python tests", "Golden tests with shipped reference data"])
def test_each_platform_runs(workflow: dict, name: str) -> None:
    assert "if" not in _step(workflow["jobs"]["test"], name)


def test_a_golden_whose_data_ship_cannot_skip(workflow: dict) -> None:
    step = _step(workflow["jobs"]["test"], "Golden tests with shipped reference data")
    assert step["env"]["POPMIX_REQUIRE_GOLDEN"] == "1"


def test_one_platform_is_enough_for_the_lints(workflow: dict) -> None:
    assert _runners(workflow["jobs"]["lint"]) == ["ubuntu-latest"]
