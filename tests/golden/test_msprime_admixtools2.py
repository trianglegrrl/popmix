"""Golden master tests: msprime scripts vs admixtools2 msprime_sim / msprime_genome.

msprime.json (generate_msprime.R) holds admixtools2's script for three
dg_6pop graphs, each with its defaults and other options, for both
functions; popmix must write the same text (its output prefix aside).
"""

from __future__ import annotations

import json
from pathlib import Path

import popmix as pm
import pytest

REF = Path(__file__).parent / "admixtools2" / "v2.0.10" / "dg_6pop" / "msprime.json"
pytestmark = pytest.mark.skipif(not REF.exists(), reason="msprime.json not generated")


def _scripts():
    if not REF.exists():
        return []
    return [pytest.param(s, id=f"{s['graph']}-{s['fun']}-{k}") for k, s in enumerate(json.loads(REF.read_text())["scripts"])]


@pytest.mark.parametrize("s", _scripts())
def test_script_text(s, tmp_path):
    fun = getattr(pm, s["fun"])
    args = s["args"] if isinstance(s["args"], dict) else {}
    outpref = tmp_path / "sim"
    path = fun([tuple(e) for e in s["edges"]], outpref=outpref, **args)
    assert Path(path).read_text() == s["text"].replace("{outpref}", str(outpref)) + "\n"


def test_r_number_format():
    from popmix.simulate import _r
    cases = {0: "0", 7000: "7000", 1e5: "1e+05", 1e6: "1e+06", 1000001: "1000001", 0.001: "0.001", 1e-4: "1e-04",
             1.25e-8: "1.25e-08", 0.7: "0.7", 1 - 0.3: "0.7", 999.0: "999", -2.5: "-2.5", 123456.5: "123456.5"}
    for x, want in cases.items():
        assert _r(x) == want, x
