"""
Popmix: Population genetics F-statistics and admixture analysis toolkit.

A Rust/Python reimplementation of admixtools2.
"""

from popmix._core import (
    # Threading
    init_threads,
    get_threads,
    version,
    __version__,
    # I/O functions
    read_plink,
    read_eigenstrat,
    read_packedancestrymap,
    read_vcf,
    read_genotypes,
    read_f2,
    dataset_info,
    list_populations,
    suggest_populations,
    convertf,
    merge,
    detect_format,
    # F-statistics functions
    extract_f2,
    f2_from_frequencies,
    # qpAdm/qpWave functions
    qpwave,
    qpadm,
    qpf4ratio,
    qpadm_rotate,
    qpwave_pairs,
    # qpGraph functions
    qpgraph,
    build_tree,
    fit_graph,
    find_graphs,
    all_graphs,
    graph_key,
    # DATES functions
    dates,
    # Classes
    PySample,
    PySnp,
    PyGenotypeMatrix,
    PyAlleleFrequencies,
    PyF2Blocks,
    PyRankTest,
    PyQpWaveResult,
    PyQpAdmResult,
    PyF4RatioResult,
    PyRotationModel,
    PyQpAdmRotateResult,
    PyQpWavePairResult,
    PyQpWavePairsResult,
    PyAdmixGraph,
    PyGraphFitResult,
    PyScoredGraph,
    PyLdBin,
    PyDecayFit,
    PyDatesResult,
)

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
from popmix.kinship import king, read_kinship
from popmix.pca import PcaResult, smartpca
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
GraphFitResult = PyGraphFitResult
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
    "joint_sfs",
    "Sfs",
    "PcaResult",
    # qpAdm/qpWave functions
    "qpwave",
    "qpadm",
    "qpf4ratio",
    "qpadm_rotate",
    "qpwave_pairs",
    # qpGraph functions
    "build_tree",
    "qpgraph",
    "fit_graph",
    "find_graphs",
    "all_graphs",
    "graph_key",
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
    "GraphFitResult",
    "ScoredGraph",
    "LdBin",
    "DecayFit",
    "DatesResult",
]
