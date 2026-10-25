"""find_graphs: the stepwise and exhaustive searches (Rust).

Their agreement with admixtools2's graphs and scores is in
tests/golden/test_find_graphs_admixtools2.py.
"""

import pytest

import popmix as pm


def test_an_unknown_method_is_an_error():
    with pytest.raises(ValueError, match="stepwise, exhaustive"):
        pm.find_graphs(None, ["A", "B", "C"], method="evolutionary")
