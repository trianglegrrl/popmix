"""read_genotypes(populations=...) reads only those samples (TASKS.md F2)."""

from pathlib import Path

import numpy as np
import pytest

import popmix as pm

FAMILY = str(Path(__file__).parents[1] / "golden" / "dates" / "family_chr20_22")


@pytest.fixture(scope="module", params=["packedancestrymap", "plink", "eigenstrat"])
def prefix(request, tmp_path_factory):
    if request.param == "packedancestrymap":
        return FAMILY
    out = str(tmp_path_factory.mktemp(request.param) / "family")
    pm.convertf(FAMILY, out, request.param)
    return out


def test_only_the_populations_samples_are_read(prefix):
    full = pm.read_genotypes(prefix)
    some = pm.read_genotypes(prefix, populations=["Yoruba", "Simulation"])
    keep = [i for i, s in enumerate(full.samples) if s.population in ("Yoruba", "Simulation")]
    assert [s.id for s in some.samples] == [full.samples[i].id for i in keep]
    assert [s.ploidy for s in some.samples] == [full.samples[i].ploidy for i in keep]
    np.testing.assert_array_equal(some.to_numpy(), full.to_numpy()[:, keep])
    assert some.n_snps == full.n_snps


def test_an_unknown_population_is_an_error():
    with pytest.raises(ValueError, match="Atlantis"):
        pm.read_genotypes(FAMILY, populations=["Atlantis"])
