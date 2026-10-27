"""No input through the Python API panics (B6).

Every call here gets hostile arguments: unknown, empty, repeated and
non-ASCII population names, NaN, infinite, negative and huge numbers, and
indices out of range. Each may raise an ordinary exception (ValueError,
TypeError, OverflowError, ...); none may raise pyo3's PanicException, which
derives from BaseException and so escapes ``except Exception``. An abort
(an allocation the size of a huge argument, say) would end the test run.
"""

from __future__ import annotations

import contextlib
import math
import os
import subprocess
import sys

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation", "Simulation.1"]
FIXTURE = "tests/fixtures/small_eigenstrat"

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")

names = st.one_of(
    st.sampled_from(POPS),
    st.sampled_from(["", "French ", "fRENCH", "ééé", "Pöp", "\0", "x" * 300, "Simulation.1.2"]),
    st.text(max_size=12),
)
name_lists = st.lists(names, max_size=5)
numbers = st.one_of(
    st.sampled_from([0.0, -0.0, -1.0, 1e-300, 1e-12, 1e300, math.nan, math.inf, -math.inf]),
    st.floats(allow_nan=True, allow_infinity=True),
)
counts = st.one_of(st.sampled_from([0, 1, 2, 255, 256, 257, 2**31, 2**63 - 1]), st.integers(-10, 10**6))

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
FEW = settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])


def survive(f, *args, **kwargs) -> None:
    """Call f; ordinary exceptions are fine, a panic is not."""
    # a PanicException derives from BaseException, so it is not suppressed
    with contextlib.suppress(Exception):
        f(*args, **kwargs)


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(SMALL, POPS)


@SETTINGS
@given(a=names, b=names, c=names, d=names)
def test_f2_statistics(f2, a, b, c, d):
    survive(f2.f2, a, b)
    survive(f2.fst, a, b)
    survive(f2.f3, a, b, c)
    survive(f2.f4, a, b, c, d)
    survive(f2.f2_per_block, a, b)
    survive(f2.subset, [a, b, c])
    survive(pm.suggest_populations, SMALL, a, limit=5)


@SETTINGS
@given(target=names, sources=name_lists, right=name_lists, fudge=numbers)
def test_qpadm_and_qpwave(f2, target, sources, right, fudge):
    survive(pm.qpadm, f2, target, sources, right, fudge=fudge)
    survive(pm.qpwave, f2, sources, right, fudge=fudge)
    survive(pm.qpadm_p, f2, sources, right, target)
    survive(pm.qpwave_pairs, f2, sources, right, fudge=fudge)
    survive(pm.qpadm_rotate, f2, sources, target)


@SETTINGS
@given(pops=st.lists(names, min_size=5, max_size=5))
def test_ratios_and_dstats(f2, pops):
    survive(pm.qpf4ratio, f2, *pops)
    survive(pm.qpdstat, f2, *pops[:4])
    survive(pm.qp3pop, f2, *pops[:3])


@FEW
@given(edges=st.lists(st.tuples(names, names), max_size=6), diag=numbers, numstart=counts)
def test_qpgraph(f2, edges, diag, numstart):
    survive(pm.qpgraph, f2, edges, diag=diag, numstart=min(numstart, 3))
    survive(pm.find_graphs, f2, [e[0] for e in edges], nadmix=min(max(numstart, 0), 1))


@SETTINGS
@given(start=counts, end=counts)
def test_snp_ranges(start, end):
    survive(pm.read_eigenstrat, FIXTURE, (start, end))
    survive(pm.read_genotypes, SMALL, (start, end))


@FEW
@given(block_size=numbers, maxmiss=numbers, pops=name_lists)
def test_extract_f2(block_size, maxmiss, pops):
    survive(pm.extract_f2, SMALL, pops or POPS[:2], block_size=block_size, maxmiss=maxmiss)


@FEW
@given(binsize=numbers, maxdis=numbers, admixed=names)
def test_ld_methods(binsize, maxdis, admixed):
    survive(pm.rolloff, SMALL, admixed, ["French", "Yoruba"], binsize=binsize, maxdis=maxdis)


@FEW
@given(k=counts, cv=counts, tolerance=numbers)
def test_admixture(k, cv, tolerance):
    survive(pm.admixture, FIXTURE, min(k, 5), cv=min(cv, 400), tolerance=tolerance, max_iterations=5)


@SETTINGS
@given(threads=counts)
def test_threads(threads):
    survive(pm.init_threads, threads)


def test_a_huge_thread_count_before_the_pool_exists():
    """rayon allocates for every thread asked for: in a fresh interpreter
    init_threads(2**40) used to abort it."""
    code = (
        "import popmix as pm\n"
        "try:\n"
        "    pm.init_threads(2**40)\n"
        "except Exception as e:\n"
        "    print(type(e).__name__, e)\n"
    )
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stderr
    assert "thread count" in run.stdout
