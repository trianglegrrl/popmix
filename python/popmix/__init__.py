"""
Popmix: Population genetics F-statistics and admixture analysis toolkit.

A Rust/Python reimplementation of admixtools2.
"""

import importlib
import sys
import types
from typing import TYPE_CHECKING

from popmix._core import (
    AdmixGraph,
    AlleleFrequencies,
    F2Blocks,
    GenotypeMatrix,
    Sample,
    Snp,
    __version__,
    all_graphs,
    build_tree,
    convertf,
    dataset_info,
    # DATES functions
    detect_format,
    # F-statistics functions
    extract_f2,
    f2_from_frequencies,
    get_threads,
    graph_key,
    # Threading
    init_threads,
    list_populations,
    merge,
    # qpGraph functions
    # qpAdm/qpWave functions
    read_eigenstrat,
    read_f2,
    read_genotypes,
    read_packedancestrymap,
    # I/O functions
    read_plink,
    read_vcf,
    suggest_populations,
    version,
)
from popmix.errors import (
    ComputationError,
    FileError,
    InputError,
    InternalError,
    MissingFileError,
    PopmixError,
    PopmixWarning,
    PopulationNotFoundError,
)

# Python submodules load on first use: most import pandas, which takes
# longer to import than the rest of popmix, and `popmix --help` needs none.
_LAZY = {
    "AdmixtureResult": "admixture",
    "DatesResult": "results",
    "DecayFit": "results",
    "F4RatioResult": "results",
    "LdBin": "results",
    "PcaResult": "pca",
    "PopDropRow": "results",
    "QpAdmAssessment": "qpadm_assess",
    "QpAdmResult": "results",
    "QpGraphEdge": "results",
    "QpGraphF3": "results",
    "QpGraphResult": "results",
    "QpGraphStart": "results",
    "QpWaveResult": "results",
    "RankTest": "results",
    "RolloffResult": "rolloff",
    "ScoredGraph": "results",
    "Sfs": "sfs",
    "admixture": "admixture",
    "admixture_loglik": "admixture",
    "boo_list": "resample",
    "cluster_order": "dstats",
    "compare_fits": "resample",
    "dates": "ld_dating",
    "dstat_matrix": "dstats",
    "extract_f2_subset": "f2_dirs",
    "f3blockdat_from_geno": "dstats",
    "f4blockdat_from_geno": "dstats",
    "find_graphs": "graph_search",
    "fstat_combinations": "dstats",
    "joint_sfs": "sfs",
    "king": "kinship",
    "loo_list": "resample",
    "outgroup_f3_matrix": "dstats",
    "parse_dot": "dot",
    "qp3pop": "dstats",
    "qp3pop_resample_inds": "resample",
    "qp3pop_resample_snps": "resample",
    "qpadm": "fits",
    "qpadm_models": "qpadm_models",
    "qpadm_multi": "qpadm_models",
    "qpadm_p": "qpadm_models",
    "qpadm_resample_inds": "resample",
    "qpadm_resample_snps": "resample",
    "qpadm_rotate": "qpadm_models",
    "qpdstat": "dstats",
    "qpdstat_resample_inds": "resample",
    "qpdstat_resample_snps": "resample",
    "qpf4ratio": "fits",
    "qpgraph": "fits",
    "qpgraph_out_of_sample": "graph_assess",
    "qpgraph_resample_inds": "resample",
    "qpgraph_resample_multi": "resample",
    "qpgraph_resample_snps": "resample",
    "qpwave": "fits",
    "qpwave_pairs": "qpadm_models",
    "qpwave_resample_inds": "resample",
    "qpwave_resample_snps": "resample",
    "read_admixtools2_f2": "f2_dirs",
    "read_kinship": "kinship",
    "rolloff": "rolloff",
    "rotate_models": "qpadm_models",
    "sample_ploidy": "dstats",
    "sfs": "sfs",
    "smartpca": "pca",
    "test_cladality": "qpadm_models",
    "write_admixtools2_f2": "f2_dirs",
    "write_dot": "dot",
}


def __getattr__(name: str) -> object:
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module 'popmix' has no attribute {name!r}")
    value = getattr(importlib.import_module(f"popmix.{module}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_LAZY))


class _Package(types.ModuleType):
    """``popmix.admixture`` is the function even once its submodule
    ``popmix.admixture`` is imported, which sets the attribute to the module."""

    def __getattribute__(self, name: str) -> object:
        if _LAZY.get(name) == name:
            return getattr(importlib.import_module(f"popmix.{name}"), name)
        return super().__getattribute__(name)


sys.modules[__name__].__class__ = _Package


if TYPE_CHECKING:
    from popmix.admixture import AdmixtureResult, admixture, admixture_loglik
    from popmix.dot import parse_dot, write_dot
    from popmix.dstats import (
        cluster_order,
        dstat_matrix,
        f3blockdat_from_geno,
        f4blockdat_from_geno,
        fstat_combinations,
        outgroup_f3_matrix,
        qp3pop,
        qpdstat,
        sample_ploidy,
    )
    from popmix.f2_dirs import extract_f2_subset, read_admixtools2_f2, write_admixtools2_f2
    from popmix.fits import qpadm, qpf4ratio, qpgraph, qpwave
    from popmix.graph_assess import qpgraph_out_of_sample
    from popmix.graph_search import find_graphs
    from popmix.kinship import king, read_kinship
    from popmix.ld_dating import dates
    from popmix.pca import PcaResult, smartpca
    from popmix.qpadm_assess import QpAdmAssessment
    from popmix.qpadm_models import (
        qpadm_models,
        qpadm_multi,
        qpadm_p,
        qpadm_rotate,
        qpwave_pairs,
        rotate_models,
        test_cladality,
    )
    from popmix.resample import (
        boo_list,
        compare_fits,
        loo_list,
        qp3pop_resample_inds,
        qp3pop_resample_snps,
        qpadm_resample_inds,
        qpadm_resample_snps,
        qpdstat_resample_inds,
        qpdstat_resample_snps,
        qpgraph_resample_inds,
        qpgraph_resample_multi,
        qpgraph_resample_snps,
        qpwave_resample_inds,
        qpwave_resample_snps,
    )
    from popmix.results import (
        DatesResult,
        DecayFit,
        F4RatioResult,
        LdBin,
        PopDropRow,
        QpAdmResult,
        QpGraphEdge,
        QpGraphF3,
        QpGraphResult,
        QpGraphStart,
        QpWaveResult,
        RankTest,
        ScoredGraph,
    )
    from popmix.rolloff import RolloffResult, rolloff
    from popmix.sfs import Sfs, joint_sfs, sfs

__all__ = [
    # Errors
    "PopmixError",
    "InputError",
    "PopulationNotFoundError",
    "MissingFileError",
    "FileError",
    "ComputationError",
    "InternalError",
    "PopmixWarning",
    # Threading
    "init_threads",
    "get_threads",
    "version",
    "__version__",
    # I/O functions
    "read_plink",
    "read_eigenstrat",
    "read_packedancestrymap",
    "read_vcf",
    "read_genotypes",
    "read_f2",
    "dataset_info",
    "list_populations",
    "suggest_populations",
    "convertf",
    "merge",
    "detect_format",
    # F-statistics functions
    "extract_f2",
    "f2_from_frequencies",
    "read_admixtools2_f2",
    "write_admixtools2_f2",
    "extract_f2_subset",
    "qpdstat",
    "qp3pop",
    "fstat_combinations",
    "f4blockdat_from_geno",
    "f3blockdat_from_geno",
    "sample_ploidy",
    "outgroup_f3_matrix",
    "dstat_matrix",
    "cluster_order",
    "smartpca",
    "king",
    "read_kinship",
    "sfs",
    "admixture",
    "rolloff",
    "RolloffResult",
    "admixture_loglik",
    "AdmixtureResult",
    "joint_sfs",
    "Sfs",
    "PcaResult",
    # qpAdm/qpWave functions
    "qpwave",
    "qpadm",
    "qpf4ratio",
    "qpadm_rotate",
    "qpwave_pairs",
    "qpadm_p",
    "qpadm_multi",
    "qpadm_models",
    "rotate_models",
    "test_cladality",
    # Resampling
    "loo_list",
    "boo_list",
    "qpgraph_resample_snps",
    "qpgraph_resample_multi",
    "compare_fits",
    "qpadm_resample_snps",
    "qpwave_resample_snps",
    "qpdstat_resample_snps",
    "qp3pop_resample_snps",
    "qpadm_resample_inds",
    "qpwave_resample_inds",
    "qpgraph_resample_inds",
    "qpdstat_resample_inds",
    "qp3pop_resample_inds",
    # qpGraph functions
    "build_tree",
    "qpgraph",
    "qpgraph_out_of_sample",
    "find_graphs",
    "all_graphs",
    "graph_key",
    "write_dot",
    "parse_dot",
    # DATES functions
    "dates",
    # Classes
    "Sample",
    "Snp",
    "GenotypeMatrix",
    "AlleleFrequencies",
    "F2Blocks",
    "RankTest",
    "QpWaveResult",
    "QpAdmResult",
    "QpAdmAssessment",
    "PopDropRow",
    "F4RatioResult",
    "AdmixGraph",
    "ScoredGraph",
    "LdBin",
    "DecayFit",
    "DatesResult",
    "QpGraphResult",
    "QpGraphEdge",
    "QpGraphF3",
    "QpGraphStart",
]
