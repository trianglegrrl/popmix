"""One input convention across the public API (TASKS.md J1).

- Options are keyword-only, so that 1.x can add or reorder them without
  breaking a call: every parameter with a default comes after ``*``, but
  the populations a statistic is of, which may be given in order as in
  admixtools2 (``qpdstat(f2, "A", "B", "C", "D")``).
- Genotype analyses take the dataset's prefix first; f2 analyses the f2
  blocks.
- The same role has the same name: qpAdm's models are ``left, right,
  target`` (admixtools2's order); the admixed population and its sources of
  an admixture date are ``admixed, sources``.
"""

from __future__ import annotations

import inspect

import pytest

import popmix as pm

# the populations of a statistic, given in order as in admixtools2
DATA = {"pop1", "pop2", "pop3", "pop4", "populations", "pca_pops", "project_pops", "projection",
        "target", "rightfix"}
# a lone optional argument where a keyword would read worse
POSITIONAL_OPTIONS = {
    "init_threads": {"threads"},
    "write_dot": {"path"},
}


def _public_functions():
    for name in pm.__all__:
        obj = getattr(pm, name)
        if callable(obj) and not isinstance(obj, type):
            yield name, obj


@pytest.mark.parametrize("name,func", list(_public_functions()), ids=lambda x: x if isinstance(x, str) else "")
def test_options_are_keyword_only(name, func):
    sig = inspect.signature(func)
    positional = {p.name for p in sig.parameters.values()
                  if p.default is not inspect.Parameter.empty
                  and p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)}
    assert positional - DATA <= POSITIONAL_OPTIONS.get(name, set())



@pytest.mark.parametrize("name", ["qpadm", "qpadm_p", "qpadm_resample_snps", "qpadm_resample_inds"])
def test_qpadm_models_are_left_right_target(name):
    params = list(inspect.signature(getattr(pm, name)).parameters)
    i = params.index("left")
    assert params[i:i + 3] == ["left", "right", "target"]


@pytest.mark.parametrize("name", ["dates", "rolloff"])
def test_admixture_dates_take_a_prefix_the_admixed_and_its_sources(name):
    assert list(inspect.signature(getattr(pm, name)).parameters)[:3] == ["prefix", "admixed", "sources"]


# J2: statistics are named as in admixtools2's tables everywhere: est, se, z,
# p, chisq, dof; these older names are gone from every result
OLD_NAMES = {"estimate", "std_err", "std_errors", "z_score", "p_value", "chi_sq", "df",
             "weight_std_errors", "is_feasible", "left_pops", "right_pops", "test_pop", "source1",
             "source2"}


def _result_classes():
    import popmix._core as core
    for name in dir(core):
        cls = getattr(core, name)
        if isinstance(cls, type) and name.startswith("Py"):
            yield name, cls


@pytest.mark.parametrize("name,cls", list(_result_classes()), ids=lambda x: x if isinstance(x, str) else "")
def test_results_use_admixtools2_names(name, cls):
    attributes = {a for a in dir(cls) if not a.startswith("_")}
    assert attributes & OLD_NAMES == set()


def test_admixture_dates_name_their_populations_as_their_inputs():
    r = pm.dates("tests/golden/dates/family_chr20_22", "Simulation.1", ["French", "Yoruba"], jackknife=False)
    assert (r.admixed, r.sources) == ("Simulation.1", ("French", "Yoruba"))
    assert list(r.to_ld_dataframe().columns)[:2] == ["distance_cm", "n_pairs"]
    assert "se" in r.to_ld_dataframe().columns
