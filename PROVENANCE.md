# Provenance

popmix is MIT-licensed (LICENSE); code derived from ADMIXTOOLS and EIGENSOFT
keeps their BSD-3 licence (NOTICE). This file records, module by module,
what each part was written from. It comes from an audit on 2026-10-25
(TASKS.md C1). The audit compared each module with its reference tool's
source side by side: structure, names, control flow and comments. The
modules it found derived from GPL or unlicensed sources were rewritten
clean-room, or removed (C2).

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

The golden tests also run dadi 2.4.4 (site frequency spectra), PLINK 1.9
(`--bmerge`, for the SNPs a merge must flip) and BCFtools (reading the VCF
popmix writes), as black boxes. KING-robust kinship was written from
Manichaikul et al. 2010 and checked against PLINK 2's implementation.
`CITATION.cff` cites every tool here, at the version its goldens ran, and
the paper of every method popmix implements (README.md, "Citation").

## Modules

Verdicts: **translation** (routine-by-routine correspondence with the same
names, control flow or comments), **partly** (mirrors structure, names or
quirks in places), **independent** (only the behaviour matches).

### Rewritten clean-room for 1.0 (C2)

These modules followed a GPL or unlicensed source too closely and were
replaced. For each, a specification of behaviour was written
(`docs/cleanroom/*.md`) from the papers, the tools' documentation and their
black-box output. A separate implementer who had not read the old code or
the source then wrote the new code from the specification, the tests and
further black-box runs. The old code was deleted unread, located by line
numbers only.

| Module | Was derived from | Specification |
|--------|------------------|---------------|
| `crates/popmix-core/src/qpadm/fit.rs` | admixtools2 C++ qpAdm fit (GPL-3) | `qpadm-fit.md` |
| `crates/popmix-core/src/qpgraph/qpgraph_fit.rs` (edge lengths, f3 covariance) | admixtools2 C++ and R (GPL-3) | `qpgraph-edge-lengths.md`, `jackknife.md` |
| `crates/popmix-core/src/rcompat/lbfgsb.rs` | R's L-BFGS-B 2.3 and `optim` (GPL-2+) | `lbfgsb.md`; now translated from BSD L-BFGS-B 3.0 (NOTICE), with 2.3's subspace truncation as behaviour |
| `crates/popmix-core/src/rcompat/rng.rs` | R's `RNG.c` (GPL-2+) | `small-helpers.md`; the Mersenne Twister from BSD `mt19937ar.c` (NOTICE) |
| `crates/popmix-core/src/fstats/jackknife.rs`, `geno_fstats.rs` (jackknife), `qpadm/f4ratio.rs` | admixtools2 R (GPL-3) | `jackknife.md`, `small-helpers.md` |
| `crates/popmix-core/src/qpgraph/graph_search.rs` (enumeration) | admixtools2 R (GPL-3) | `small-helpers.md`; the exclusion rule and the order were found from `generate_all_graphs` output |
| `python/popmix/resample.py` (`compare_fits`), `qpadm_models.py` (model sets), `dstats.py` (`fstat_combinations`) | admixtools2 R (GPL-3) | `small-helpers.md` |
| `crates/popmix-core/src/dates/` (the whole module: weights, LD decay, mesh, exponential fit, jackknife, analysis) | DATES v753 `dates.c` (no licence) | `dates-mesh.md`, then `dates.md` |

Exposure during the rewrites, recorded for completeness:
- The f4ratio implementer had earlier seen about 25 lines of the old body while reading a caller of the jackknife.
- The graph-enumeration implementer saw the old three-line power-set helper and two closing lines of the old enumeration in range prints.
- One R traceback from a black-box run showed internal names.

None of these is the substance of the new code.

### BSD-3: kept, with attribution (NOTICE)

| Module | From |
|--------|------|
| `crates/popmix-core/src/rolloff.rs` | ADMIXTOOLS `rolloff.c`, `ldsubs.c` |
| `crates/popmix-core/src/pca/tracy_widom.rs`, `tw_table.rs` | EIGENSOFT `twstats`, `twtable.h` |
| `crates/popmix-core/src/pca/smartpca.rs` (partly: names) | EIGENSOFT `smartpca.c` |

### Independent

Everything else: the f-statistics, f2 blocks and blocks, qpWave, qpAdm,
qpadm_p, rotation and pairs, SFS, genotype matrices, KING-robust kinship
(which differs in implementation from PLINK 2's), ADMIXTURE
(SQUAREM, where ADMIXTURE's source is unpublished), every reader and writer in `popmix-io` (merge, VCF, TGENO),
and the Python bindings, CLI, f2 directories, DOT, RDS, plots and SFS
wrappers.

### Corrected by the independent review (2026-10-29, C7)

- **READ kinship** (`crates/popmix-core/src/kinship/read.rs`): **partly**.
  It was written (0733e04, before the provenance rule) with READv2's GPL-3
  source in view: it named READv2's functions, and DEVIATIONS.md quoted its
  code. Rewritten clean-room from `docs/cleanroom/read.md` (TASKS.md C9).
- **qpGraph's admixture-weight fit and path weights**
  (`qpgraph/qpgraph_fit.rs`, d766673 and 473ed7d): written before the rule
  and given no verdict by the C1 audit, which covered the edge lengths and
  the f3 covariance (both rewritten in C2). Its comments now describe
  behaviour; it is left to the similarity review below.
- **Compat flags first found in admixtools2's source.** ADMIXTOOLS2-0001
  (popdrop's path), -0003 (`joint_sfs` counts), -0004 (`write_dot` names) and
  -0006 (chromosome codes) were found by reading admixtools2's R code
  before the rule. The flags reproduce observable output and no code was
  copied; their DEVIATIONS.md entries now describe the behaviour, and every
  later entry was found from black-box runs. The same holds for the RNG's
  edge clamp in `docs/cleanroom/small-helpers.md`.
- **`python/popmix/rds.py`**: written from the R Internals manual's
  description of serialization and checked against admixtools2's f2
  directories; not compared with R's source.

### How popmix was written

popmix was written by its author with AI coding agents (Anthropic's Claude,
through Claude Code); most commits are co-authored by one. The clean-room
rewrites (C2, C9) followed the process above: one agent wrote the
specification from papers, documentation and black-box runs, and another,
which had not seen the old code or the reference's source, implemented it.
Before the provenance rule (2026-10-25), agents translated some modules
from the reference tools' source: those are the modules rewritten in C2 or
removed (above). Since the rule, no agent works from a GPL tool's source.
The models were
trained on public code, which may include the reference tools; a
similarity comparison of popmix against those sources, by someone outside
development, is a question for the legal review (C7).

## Removed

- **ALDER** (`alder/`, the `alder` function and command): a translation of
  ALDER 1.03, which is licensed for academic and non-profit use only.
  Removed for 1.0; DATES and ROLLOFF date admixture.
- **Translations of admixtools2 (GPL-3), removed for 1.0** and to return as
  rewrites: the evolutionary `find_graphs` search and its moves (the
  stepwise and exhaustive searches, popmix's own, remain) and the graph
  constraints (`satisfies_*`, `summarize_*`); `lazadm`; `qpfstats` and
  `qpf4diff`; the msprime script generators and `pseudo_dates`; the parsers
  of classic ADMIXTOOLS output (`parse_qp*`). The graph plots' layout was
  rewritten (it used `pseudo_dates`).

## Test scripts under GPL-3

`tests/golden/admixtools2/generate_reference.R`, `generate_all_graphs.R`,
`generate_dstat.R` and `generate_sfs.R` run admixtools2 to make the golden
references, and call its internal functions (`admixtools:::`) for values
its public API does not return. They are GPL-3.0-or-later like admixtools2
(SPDX headers, `tests/golden/admixtools2/COPYING`). They are test tooling,
in no popmix package; `tests/python/test_licence_marks.py` keeps it so.

## Test patches

`tests/golden/{dates,rolloff,eigensoft}/*.patch` modify the reference
tools to print full precision for the golden tests, and carry a few lines of
their source as context. They are test tooling, not part of the packages
(neither the sdist nor the wheels include `tests/`).
