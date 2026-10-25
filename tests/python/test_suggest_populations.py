"""suggest_populations() takes any text (a non-ASCII query used to panic)."""

import popmix as pm

FIXTURE = "tests/fixtures/small_eigenstrat"


def test_non_ascii_query():
    assert isinstance(pm.suggest_populations(FIXTURE, "ééé", 5), list)
    assert isinstance(pm.suggest_populations(FIXTURE, "Pöp1x", 5), list)


def test_three_character_match_counts_characters():
    names = [p for p, _ in pm.suggest_populations(FIXTURE, "popzzzz", 5)]
    assert names and all(n.lower().startswith("pop") for n in names[:1])
