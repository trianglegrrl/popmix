"""Tests for popmix.rds (R serialization). Reading R's own files is checked
in tests/golden/test_f2_dirs_admixtools2.py, on admixtools2's f2 directory."""

import math

import pytest

from popmix.rds import RObject, matrix, read_rds, write_rds


@pytest.mark.parametrize("obj", [
    [1.5, -2.0, 1e-300],
    [582, 976, None],
    ["Han.DG", None, "é"],
    [True, False],
    RObject([["x", None], [1.0, 2.0]], {"names": ["a", "b"]}),
    None,
])
def test_round_trip(tmp_path, obj):
    f = tmp_path / "x.rds"
    write_rds(str(f), obj)
    assert read_rds(str(f)) == obj


def test_matrix_and_nan(tmp_path):
    f = tmp_path / "m.rds"
    write_rds(str(f), matrix([[0.25, math.nan], [5.0, 6.0]], ["f2", "counts"]), compress=False)
    m = read_rds(str(f))
    assert m.attributes == {"dim": [2, 2], "dimnames": [None, ["f2", "counts"]]}
    assert m.data[0] == 0.25 and math.isnan(m.data[1]) and m.data[2:] == [5.0, 6.0]


def test_not_rds(tmp_path):
    f = tmp_path / "bad.rds"
    f.write_bytes(b"hello")
    with pytest.raises(ValueError, match="not an XDR"):
        read_rds(str(f))
