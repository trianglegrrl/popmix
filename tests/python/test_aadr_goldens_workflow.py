"""The weekly AADR golden workflow fetches the Zenodo subset and runs the
goldens on it (TASKS.md K2)."""

from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/aadr-goldens.yml"


@pytest.fixture(scope="module")
def workflow() -> dict:
    return yaml.safe_load(WORKFLOW.read_text())


@pytest.fixture(scope="module")
def steps(workflow: dict) -> list[dict]:
    return workflow["jobs"]["goldens"]["steps"]


def _run(steps: list[dict]) -> str:
    return "\n".join(s.get("run", "") for s in steps)


def test_runs_weekly_and_by_hand_only(workflow: dict) -> None:
    on = workflow.get("on", workflow.get(True))  # YAML 1.1 reads `on` as True
    assert set(on) == {"schedule", "workflow_dispatch"}


def test_the_schedule_is_once_a_week(workflow: dict) -> None:
    on = workflow.get("on", workflow.get(True))
    minute, hour, day, month, weekday = on["schedule"][0]["cron"].split()
    assert (day, month, weekday.isdigit()) == ("*", "*", True)


def test_fetches_the_subset_with_the_archive_script(steps: list[dict]) -> None:
    assert "tests/golden/aadr/archive_aadr_v66p1.py fetch" in _run(steps)


def test_caches_the_subset_keyed_on_the_archive_script(steps: list[dict]) -> None:
    cache = [s for s in steps if s.get("uses", "").startswith("actions/cache")]
    assert "hashFiles('tests/golden/aadr/archive_aadr_v66p1.py')" in cache[0]["with"]["key"]


def test_runs_the_golden_tests(steps: list[dict]) -> None:
    assert "pytest tests/golden --run-slow" in _run(steps)


def test_reads_the_subset_where_the_goldens_look(workflow: dict) -> None:
    assert workflow["env"]["POPMIX_REFERENCE_CACHE"] == "${{ github.workspace }}/.references"
