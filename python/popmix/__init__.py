"""
Popmix: Population genetics F-statistics and admixture analysis toolkit.

A Rust/Python reimplementation of admixtools2.
"""

from popmix._core import (
    PyAdmixGraph,
    PyAlleleFrequencies,
    PyDatesResult,
    PyDecayFit,
    PyF2Blocks,
    PyF4RatioResult,
    PyGenotypeMatrix,
    PyLdBin,
    PyQpAdmResult,
    PyQpAdmRotateResult,
    PyQpWavePairResult,
    PyQpWavePairsResult,
    PyQpWaveResult,
    PyRankTest,
    PyRotationModel,
    # Classes
    PySample,
    PyScoredGraph,
    PySnp,
    __version__,
    all_graphs,
    build_tree,
    convertf,
    dataset_info,
    # DATES functions
    dates,
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
    qpadm,
    qpf4ratio,
    # qpGraph functions
    qpgraph,
    # qpAdm/qpWave functions
    qpwave,
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
from popmix.graph_search import find_graphs
from popmix.kinship import king, read_kinship
from popmix.pca import PcaResult, smartpca
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
from popmix.rolloff import RolloffResult, rolloff
from popmix.sfs import Sfs, joint_sfs, sfs

# Convenient aliases
Sample = PySample
Snp = PySnp
GenotypeMatrix = PyGenotypeMatrix
AlleleFrequencies = PyAlleleFrequencies
F2Blocks = PyF2Blocks
RankTest = PyRankTest
QpWaveResult = PyQpWaveResult
QpAdmResult = PyQpAdmResult
F4RatioResult = PyF4RatioResult
RotationModel = PyRotationModel
QpAdmRotateResult = PyQpAdmRotateResult
QpWavePairResult = PyQpWavePairResult
QpWavePairsResult = PyQpWavePairsResult
AdmixGraph = PyAdmixGraph
ScoredGraph = PyScoredGraph
LdBin = PyLdBin
DecayFit = PyDecayFit
DatesResult = PyDatesResult

__all__ = [
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
    "F4RatioResult",
    "RotationModel",
    "QpAdmRotateResult",
    "QpWavePairResult",
    "QpWavePairsResult",
    "AdmixGraph",
    "ScoredGraph",
    "LdBin",
    "DecayFit",
    "DatesResult",
]
