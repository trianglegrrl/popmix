"""The classes popmix returns are named, immutable and picklable (TASKS.md J4).

Every public class is popmix's by module and has its public name (no
``builtins.PyQpAdmResult`` in a repr). Results of analyses are frozen
dataclasses; the data containers (F2Blocks, GenotypeMatrix, ...) are
Rust classes that cannot be changed in place. All of them pickle, so they
can go to and from worker processes.
"""

from __future__ import annotations

import dataclasses
import pickle

import numpy as np
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"
POPS = ["French", "Yoruba", "Simulation.1", "Simulation"]


def _public_classes():
    for name in pm.__all__:
        obj = getattr(pm, name)
        if isinstance(obj, type) and not issubclass(obj, (BaseException, Warning)):
            yield name, obj


@pytest.mark.parametrize("name,cls", list(_public_classes()), ids=lambda x: x if isinstance(x, str) else "")
def test_public_classes_are_popmix_s_under_their_public_names(name, cls):
    assert cls.__name__ == name
    assert cls.__module__.startswith("popmix")


@pytest.fixture(scope="module")
def f2():
    return pm.extract_f2(SMALL, POPS)


@pytest.fixture(scope="module")
def results(f2):
    graph = [("R", "Simulation.1"), ("R", "n1"), ("n1", "French"), ("n1", "n2"), ("n2", "Yoruba"),
             ("n2", "Simulation")]
    return {
        "qpwave": pm.qpwave(f2, ["French", "Yoruba"], ["Simulation", "Simulation.1"]),
        "qpadm": pm.qpadm(f2, ["French"], ["Yoruba", "Simulation.1"], "Simulation"),
        "qpf4ratio": pm.qpf4ratio(f2, "Simulation", "French", "Yoruba", "Simulation.1", "French"),
        "qpgraph": pm.qpgraph(f2, graph, numstart=2),
        "dates": pm.dates(SMALL, "Simulation.1", ["French", "Yoruba"], jackknife=False),
    }


@pytest.mark.parametrize("name", ["qpwave", "qpadm", "qpf4ratio", "qpgraph", "dates"])
def test_analysis_results_are_frozen_dataclasses_that_pickle(results, name):
    r = results[name]
    assert dataclasses.is_dataclass(r)
    assert type(r).__module__ == "popmix.results"
    field = dataclasses.fields(r)[0].name
    with pytest.raises(dataclasses.FrozenInstanceError):
        setattr(r, field, None)
    back = pickle.loads(pickle.dumps(r))
    assert _same(back, r)
    assert not repr(r).startswith("<builtins")


def _same(a, b) -> bool:
    if dataclasses.is_dataclass(a):
        return type(a) is type(b) and all(_same(getattr(a, f.name), getattr(b, f.name))
                                          for f in dataclasses.fields(a))
    if isinstance(a, np.ndarray):
        return np.array_equal(a, b, equal_nan=a.dtype.kind == "f")
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, float) and np.isnan(a):
        return np.isnan(b)
    return a == b


def test_f2_blocks_pickle_to_the_same_statistics(f2):
    back = pickle.loads(pickle.dumps(f2))
    assert type(back) is pm.F2Blocks
    assert back.populations == f2.populations
    assert back.block_lengths == f2.block_lengths
    assert back.f3("Simulation", "French", "Yoruba") == f2.f3("Simulation", "French", "Yoruba")
    assert back.fst("French", "Yoruba") == f2.fst("French", "Yoruba")


def test_genotypes_pickle(f2):
    geno = pm.read_genotypes(SMALL, populations=["French", "Yoruba"])
    back = pickle.loads(pickle.dumps(geno))
    assert np.array_equal(back.to_numpy(), geno.to_numpy())
    assert [(s.id, s.population, s.ploidy) for s in back.samples] == \
        [(s.id, s.population, s.ploidy) for s in geno.samples]
    assert [(s.id, s.chromosome, s.position_bp, s.position_cm) for s in back.snps[:50]] == \
        [(s.id, s.chromosome, s.position_bp, s.position_cm) for s in geno.snps[:50]]
    sample = pickle.loads(pickle.dumps(geno.samples[0]))
    assert (sample.id, sample.population) == (geno.samples[0].id, geno.samples[0].population)


def test_data_containers_cannot_be_changed_in_place(f2):
    with pytest.raises(AttributeError):
        f2.populations = []
