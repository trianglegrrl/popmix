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
