"""Tests for the qpadm() compat argument (DEVIATIONS.md, ADMIXTOOLS2-0001)."""

import os

import pytest

import popmix as pm

PREFIX = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(
    not os.path.exists(f"{PREFIX}.geno"), reason="test data not available"
)


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(PREFIX, populations=["Simulation", "French", "Yoruba", "Simulation.1"],
                         auto_only=False)


def test_unknown_compat_is_rejected(f2):
    with pytest.raises(ValueError, match="admixtools2-2.0.10"):
        pm.qpadm(f2, ["French"], ["Yoruba", "Simulation.1"], "Simulation", compat="admixtools2-9")


def test_compat_leaves_a_one_source_model_unchanged(f2):
    """With one source there is no nested path, so only `best` on the full
    model can differ (popmix puts it on the path; admixtools2 marks nothing)."""
    args = (f2, ["French"], ["Yoruba", "Simulation.1"], "Simulation")
    fixed, compat = pm.qpadm(*args), pm.qpadm(*args, compat="admixtools2-2.0.10")
    assert fixed.weights == compat.weights
    assert fixed.chisq == compat.chisq
    assert [(x.pattern, x.best) for x in fixed.popdrop] == [("0", True)]
    assert [(x.pattern, x.best) for x in compat.popdrop] == [("0", False)]



LEFT, RIGHT, TARGET = ["French_1", "Yoruba_1", "Simulation"], ["French_2", "Yoruba_2", "French_3", "Yoruba_3"], \
    "Simulation.1"


@pytest.fixture(scope="module")
def family6(tmp_path_factory) -> str:
    """The skill examples' family6: eight populations, so three sources fit."""
    from test_skill import example_data

    d = tmp_path_factory.mktemp("family6")
    example_data(d)
    return str(d / "family6")


def _cli(prefix: str, *extra: str) -> str:
    import subprocess

    return subprocess.run(["popmix", "qpadm", "-i", prefix, "--target", TARGET, "--left", ",".join(LEFT),
                           "--right", ",".join(RIGHT), *extra], capture_output=True, text=True, check=True).stdout


@pytest.mark.parametrize("compat", [None, "admixtools2-2.0.10"])
def test_the_cli_gives_the_popdrop_table_and_takes_compat(family6, compat):
    """The CLI's qpadm JSON has the popdrop table, and --compat reproduces
    admixtools2's (DEVIATIONS.md step 2: every compat flag on the CLI)."""
    import json

    r = pm.qpadm(pm.extract_f2(family6, LEFT + RIGHT + [TARGET]), LEFT, RIGHT, TARGET, compat=compat)
    rows = json.loads(_cli(family6, "--format", "json", *(["--compat", compat] if compat else [])))["popdrop"]
    assert [row["pattern"] for row in rows] == [x.pattern for x in r.popdrop]
    assert [row["best"] for row in rows] == [x.best for x in r.popdrop]
    assert [row["p"] for row in rows] == pytest.approx([x.p for x in r.popdrop], rel=1e-12)


def test_compat_changes_the_cli_s_nested_path(family6):
    import json

    fixed = json.loads(_cli(family6, "--format", "json"))["popdrop"]
    compat = json.loads(_cli(family6, "--format", "json", "--compat", "admixtools2-2.0.10"))["popdrop"]
    assert [row["best"] for row in fixed] != [row["best"] for row in compat]


def test_the_cli_text_shows_the_popdrop_table(family6):
    out = _cli(family6)
    r = pm.qpadm(pm.extract_f2(family6, LEFT + RIGHT + [TARGET]), LEFT, RIGHT, TARGET)
    assert "pattern" in out
    for row in r.popdrop:
        assert f"\n  {row.pattern} " in out
