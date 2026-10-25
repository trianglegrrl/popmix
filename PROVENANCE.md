# Provenance

popmix is MIT-licensed (LICENSE); code derived from ADMIXTOOLS and EIGENSOFT
keeps their BSD-3 licence (NOTICE). This file records, module by module,
what each part was written from, and what must change before a release
under those licences. It comes from an audit on 2026-10-25 (TASKS.md C1)
that compared each module with its reference tool's source: structure,
names, control flow and comments, side by side.

Matching a reference tool's numbers exactly is the point of popmix and needs
no permission: algorithms, formulas, file formats and outputs are not
protected. Translating a tool's source is: a translation is a derivative
work, under the tool's licence. Since the audit, development works from
papers, documentation and black-box output only (CLAUDE.md, "Provenance").

## Reference tools and their licences

| Tool | Licence | Where recorded |
|------|---------|----------------|
| admixtools2 2.0.10 | GPL-3 | its DESCRIPTION |
| R (`src/appl/lbfgsb.c`, `stats/src/optim.c`, `src/main/RNG.c`) | GPL-2-or-later | file headers |
| L-BFGS-B 3.0 (Zhu, Byrd, Lu, Nocedal; Morales, Nocedal) | BSD-3 | its distribution's README |
| MT19937 reference code (Matsumoto, Nishimura) | BSD-3 | `mt19937ar.c` header |
| ADMIXTOOLS (ROLLOFF, qp* programs, convertf, mergeit) | BSD-3, Broad Institute | `src/LICENSE.txt` |
| EIGENSOFT (smartpca, twstats) | BSD-3, Broad Institute and Harvard Medical School | `LICENSE.txt` |
| DATES v753 | none stated for its own files (its `nicksrc/` is BSD-3, Broad) | no LICENSE file or headers |
| ALDER 1.03 | academic and non-profit use only | `README.txt` §7 |
| PLINK 2 (plink-ng) | GPL-3 | `COPYING` |
| READv2 | GPL-3 | `LICENSE` |
| ADMIXTURE 1.3.0 | closed source, free for academic use | distribution |

## Modules

Verdicts: **translation** (routine-by-routine correspondence with the same
names, control flow or comments), **partly** (mirrors structure, names or
quirks in places), **independent** (only the behaviour matches).

### Must be rewritten (GPL or non-free source)

| Module | Translated from | Licence | Work |
|--------|-----------------|---------|------|
| `crates/popmix-core/src/qpadm/fit.rs` | admixtools2 C++ qpAdm fit (alternating least squares, weights) | GPL-3 | from Patterson et al. 2012, Harney et al. 2021 |
| `crates/popmix-core/src/qpadm/lazadm.rs` | admixtools2 `lazadm` | GPL-3 | from its documentation, or drop |
| `crates/popmix-core/src/qpgraph/qpgraph_fit.rs` (edge lengths; the weight-index helper partly) | admixtools2 C++ edge-length fit | GPL-3 | from Patterson et al. 2012, Maier et al. 2023 |
| `python/popmix/graph_search.py` | admixtools2 `find_graphs` and its moves | GPL-3 | from the documentation and Maier et al. 2023, own moves |
| `python/popmix/qpfstats.py` | admixtools2 `qpfstats` | GPL-3 | from Maier et al. 2023 |
| `python/popmix/graph_constraints.py` | admixtools2 graph constraint checks | GPL-3 | from the documented definitions |
| `python/popmix/simulate.py` | admixtools2 msprime script generator | GPL-3 | from the msprime documentation |
| `python/popmix/classic.py` | admixtools2 output parsers | GPL-3 | against sample ADMIXTOOLS output |
| `crates/popmix-core/src/alder/` (all four files) | ALDER 1.03 | academic only | from Loh et al. 2013 |
| `crates/popmix-core/src/rcompat/lbfgsb.rs` | R's L-BFGS-B 2.3 (f2c) and `optim` driver | GPL-2+ | re-base the core on BSD L-BFGS-B 3.0 (keeping 2.3's behaviour, as behaviour), rewrite the driver |

### Partly; to rewrite with the above

| Module | Mirrors | Licence |
|--------|---------|---------|
| `crates/popmix-core/src/qpgraph/graph_search.rs` | admixtools2's graph enumeration order (to give the same graphs in the same order) | GPL-3 |
| `crates/popmix-core/src/fstats/jackknife.rs`, `geno_fstats.rs` (helpers), `qpadm/f4ratio.rs` | admixtools2's jackknife steps, reproduced literally | GPL-3 |
| `python/popmix/resample.py`, `qpadm_models.py`, `dstats.py` (combination order) | admixtools2 helpers | GPL-3 |
| `crates/popmix-core/src/rcompat/rng.rs` | R's `RNG.c` around the Mersenne Twister | GPL-2+ (the twister is BSD MT19937) |
| `crates/popmix-core/src/dates/ld_decay.rs` (mesh tags, lag-sum binning) | DATES `dates.c` | none stated |

### BSD-3: kept, with attribution (NOTICE)

| Module | From |
|--------|------|
| `crates/popmix-core/src/rolloff.rs` | ADMIXTOOLS `rolloff.c`, `ldsubs.c` |
| `crates/popmix-core/src/pca/tracy_widom.rs`, `tw_table.rs` | EIGENSOFT `twstats`, `twtable.h` |
| `crates/popmix-core/src/pca/smartpca.rs` (partly: names) | EIGENSOFT `smartpca.c` |

### Independent

Everything else: the f-statistics, f2 blocks and blocks, qpWave, qpAdm,
qpadm_p, rotation and pairs, SFS, genotype matrices, kinship (KING-robust
and READ, which differ in implementation from PLINK 2 and READv2), ADMIXTURE
(SQUAREM, where ADMIXTURE's source is unpublished), DATES apart from the two
helpers above, every reader and writer in `popmix-io` (merge, VCF, TGENO),
and the Python bindings, CLI, f2 directories, DOT, RDS, plots and SFS
wrappers.

## Test patches

`tests/golden/{alder,dates,rolloff,eigensoft}/*.patch` modify the reference
tools to print full precision for the golden tests, and carry a few lines of
their source as context. They are test tooling, not part of the packages
(neither the sdist nor the wheels include `tests/`).
