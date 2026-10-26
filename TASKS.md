# Review tasks

Findings of the five-reviewer pre-release review of commit fc7d00c
(2026-10-25), plus the decisions that came out of it. Work through them one
at a time, TDD as always; tick a box when its fix and test are committed and
note the commit.

Reviewers: **PG** population-genetics bioinformatician, **PR** emeritus
professor of computational popgen, **PD** CS postdoc (genomics MSc),
**AI** AI-built-code expert, **RS** senior Rust programmer.

Status of each finding: **[confirmed]** reproduced again outside the review,
**[reviewer-verified]** the reviewer ran or measured it, **[suspected]** read
from the code, not yet shown.

## Decisions

- Licence stays **MIT** (2026-10-25). Section C makes that defensible.
- **Performance and scalability are requirements**, not nice-to-haves:
  popmix must be fast, flexible, and use what safe Rust can buy (section F).
- Results against the reference tools, with timings and their variation,
  must be visible: a benchmark report (section E).
- The 1.0.0 release waits for sections A to E and O.
- Reference data is **AADR v66.p1** from now on (section O). It removes 161
  modern individuals of Jacobs et al. 2019 (Cell 177:1010) that v62.0 and
  v66.0 contained; none of them, nor their 27 populations, appear in the
  repository, the reference cache or the nf-core test data (checked
  2026-10-25).
- popmix ships an agent skill documenting every CLI command (section N).
- `popmix convertf` becomes a full replacement for EIGENSOFT/AdmixTools
  convertf (section P).

## A. Wrong results (P0)

- [x] **A1. VCF genotypes inverted on conversion and merge.** *Done in 51c2f98.* `vcf.rs:113-118`
  stores the ALT dosage with column 5 = REF, but an EIGENSTRAT digit counts
  copies of column 5; `merge.rs:86-105` then "flips" on allele letters. A
  REF=A sample `0/0` becomes GG in EIGENSTRAT and PLINK output (plink2 reads
  it back as GG); merging VCF with AADR mirrors the user's samples.
  PG [confirmed]
- [x] **A2. Pseudohaploid samples not detected outside packed EIGENSTRAT.** *Done in a22fb64.*
  `popmix-python/src/lib.rs:1466-1484` never calls `detect_ploidy` for
  PLINK, text EIGENSTRAT, VCF or PED; nor does `read_genotypes()
  .allele_frequencies()` (README "Method 1"). Same samples: f2(A, Anc) 0.0206
  packed vs 0.0485 text/PLINK. Affects f2, FST, f3, qpGraph for aDNA.
  admixtools2 corrects PLINK input. PG [confirmed]
- [x] **A3. f3/f4 wrong when a population repeats.** *Done in e3d4f45.* `fstats/f3_f4.rs:75-86,
  ~132` read the self-pair block value −2p(1−p)/(n−1) (the comment in
  `fstats/f2.rs:58-69` claiming f2(A,A)=0 is false); admixtools2 zeroes the
  diagonal. f3(French; French, Yoruba) = −0.00288, z = −177
  (should be 0); f4(F,Y;F,S) 6% off f3(F;Y,S). Check qpdstat/qp3pop with
  `unique_only=False`. PD [confirmed]
- [x] **A4. Saturated qpAdm (dof 0) reported as p = 1, "Good model fit".** *Done in 982142e.*
  `qpadm/fit.rs:464-467`; admixtools2 gives p = NA (checked: French.DG =
  Han.DG + Yoruba.DG, right Mbuti.DG, Chimp.REF, weights 1.22/−0.22).
  Not in DEVIATIONS.md. PR, PG [confirmed]
- [x] **A5. Unknown, duplicated or overlapping populations accepted.** *Done in 6190214.*
  - `extract_f2(..., ["French","Yoruba","NoSuchPop"])` succeeds; every
    statistic with it returns `None` [confirmed]. Other paths drop the name or
    give `n_blocks == 0` (`genotype_matrix.rs:267`, `streaming.rs:68`)
    PD, RS, PG [reviewer-verified]
  - duplicated population overwrites the first column (`streaming.rs`
    `pop_to_col.insert`) PD [reviewer-verified]
  - qpAdm with target or a source also in right runs and reports a good fit
    (weights −235/+236); admixtools2 stops with "Duplicated pops"; the CLI de-duplicates silently (`cli/__init__.py:138`)
    PD, PG [reviewer-verified]
  - `sources=['French','French']` fails with a cryptic singular-system error
    PD [reviewer-verified]
  - `filter_populations` documented as "silently ignored"; F2Blocks
    `f2/f3/f4/fst` return `Option`, merging "unknown" with "no data" RS
- [x] **A6. Ploidy rule misclassifies diploids.** *Done in 820af03.* "Pseudohaploid if no het in
  the first 1,000 SNPs" (DEVIATIONS.md "Matched, not corrected") halves the
  hets of AADR NA18606.DG, a high-coverage genome, changing f2. By our own
  policy: fix by default (as `sfs` `pseudohaploid="detect"`), old rule under
  compat. PR [reviewer-verified from docs]
- [x] **A7. Multiallelic VCF records miscoded.** *Done in aede4ec.* `vcf.rs:110-118` sums allele
  indices: ALT=G,T with GT 0/2 (A/T het) written as homozygous G; indels kept
  as one-letter SNPs (AT/A → `A A`); ALT='.' sites kept; no warning.
  PG [reviewer-verified]
- [x] **A8. `f2_matrix()` disagrees with `f2()`.** *Done in 8493780.* `f2_blocks.rs:550-565`
  unweighted mean of block means vs the length-weighted estimate.
  PD [suspected, from code]
- [x] **A9. DATES/ALDER sums may depend on thread count.** *Done in 99266cc.* `rayon::reduce`
  over floats in `dates/ld_decay.rs:304,324`, `alder/ldcorr.rs:100` (the
  pattern fixed for PCA in 4ab0928). Identical at 1/3/4/16 threads on 4 test
  samples, too few to rule out. PD [suspected]
- [x] **A10. Merge drops or corrupts SNPs silently.** *Done in dbf87e2.* No strand flipping
  (opposite-strand SNPs dropped: 13,333 of 20,000 kept, reason not printed;
  `snps_mismatch`/`snps_strand_ambiguous` counted, never shown); PLINK `0`
  allele treated as mismatch, so a single new pseudohaploid sample gives
  "No overlapping SNPs"; duplicate (chr, pos) in file B overwrite each other
  (`merge.rs:125-128`). PG [reviewer-verified; duplicates suspected]

## B. Crashes and input validation (P0)

- [x] **B1. Out-of-range `snp_range` aborts Python.** *Done in 86cdd8e.* `popmix-io/src/plink.rs:142`
  `read_snp_range` (and `read_snp(idx)`) never checks bounds:
  `read_plink(HO_public, (0, 10**9))` → uncatchable 72 GB allocation abort;
  smaller values panic. RS [reviewer-verified]
- [x] **B2. Bad options panic.** *Done in cc6fad0.* `rolloff.rs:388` `binsize=0` → capacity
  overflow `PanicException`; tiny binsize → OOM abort. Validate options up
  front; sweep the ~15 unchecked `f64 as usize` casts (`ld_decay.rs:166,236,
  312,468,496,554,609`, `alder/mod.rs:248,264`, `f3_f4.rs:189,252`,
  `sfs.rs:146`, `streaming.rs:185`). RS [reviewer-verified]
- [x] **B3. Non-ASCII query panics.** *Done in 3513377.* `popmix-python/src/lib.rs:952`
  `&query_lower[..3]`; `suggest_populations(P, "ééé", 5)`. RS [confirmed]
- [x] **B4. NaN panics.** *Done in c5efc9c.* (`qpwave_pairs.rs:54` and the rotation sort fixed
  with A4.) `partial_cmp(..).unwrap()` at `qpwave_pairs.rs:54`,
  `qpgraph/fitting.rs:85`, `dates/ld_decay.rs:275` (sorts input positions),
  `lib.rs:2082,2348`; use `f64::total_cmp`. RS, AI [suspected] (All four
  sites fixed; no `partial_cmp(..).unwrap()` is left.)
- [x] **B5. Errors swallowed in model loops.** *Done in 8771f69.* (Names are checked before the
  loops since A5, so a typo is an error; the `Err(_) => continue` remain.) `qpadm/rotate.rs:156`,
  `qpadm/qpwave_pairs.rs:110` `Err(_) => continue`: a typo returns `Ok` with
  zero models. Validate names first; skip only a typed singular-matrix error.
  RS [reviewer-verified] (Rotation also skips, and counts, models with more
  sources than outgroups, which cannot be tested.)
- [x] **B6. A test that no user input can panic or abort through the Python
  API** *Done in 12622c4.* (fuzz/property test over the bindings). RS
- [x] **B7. `.ind`/`.fam` count not checked against the genotype file.** *Done in d5ddb03.* The
  packed `.geno` header (`GENO 28 20000`) is never read; deleting one `.ind`
  line silently relabels every later sample (f2 0.0292 → 0.0230). Same for
  PLINK `.fam`. PG [reviewer-verified] (A `.bed` has no counts: its size
  check catches a lost `.fam` line only when the record length changes.)
- [x] **B8. Text-format parsing defaults.** *Done in 2a11b67.* `eigenstrat.rs:294` unparsable cM
  becomes 0; `.snp` must have 6 columns (EIGENSOFT accepts 4); alleles
  truncated to one character. PG [reviewer-verified] (convertf 8 also
  requires 6 columns and keeps an allele's first character: matched. Only
  the genetic position is now stricter than convertf.)
- [x] **B9. Unsorted SNPs give wrong blocks.** *Done in f829c01.* `fstats/blocks.rs:77-97` never
  checks order; a decreasing position does not start a block. PD, PG
  [suspected]
- [x] **B10. VCF haploid calls and FILTER.** *Done in d8b6a5a.* Haploid `1` (male chrX) becomes
  missing; FILTER ignored. Only matters with `--all-chromosomes`.
  PG [reviewer-verified]
- [x] **B11. Graph search aborted by one singular candidate** *Done in 4487bec.*
  (`graph_search.rs` `collect::<Result<_>>`). PD [suspected, from code]
- [x] **B12. CLI flattens errors.** *Done in c4aab05.* ~25 `except Exception as e: echo; exit 1`
  and `except Exception: pass` at `cli/__init__.py:73`. AI
- [x] **B13. `admixture_cv` folds truncated to u8** *Done in 0ad353c.* (`admixture.rs:741`) when
  folds > 256. RS

## C. Licence and provenance (P0, for MIT)

Research (2026-10-25, Perplexity plus primary sources):
- A close translation of source code is an adaptation/derivative work: US
  17 U.S.C. §101 lists "translation"; EU Software Directive 2009/24 art. 4
  covers "translation, adaptation"; the GPL FAQ treats translation as
  modification. Changing language does not escape it, and AI translation
  gets no exemption.
- Not protected: algorithms, functionality, numerical output, file formats,
  languages (SAS v World Programming, CJEU C-406/10 and EWCA 2013; Altai
  abstraction-filtration-comparison with merger and scènes à faire; Google v
  Oracle was a fair-use holding about interfaces, not a licence to translate
  implementations). Matching admixtools2's numbers exactly is fine; copying
  its code structure, decomposition, names and comments is the risk.
- Practice for MIT rewrites of GPL tools (uutils/coreutils): implement from
  papers, documentation and black-box behaviour, never from the GPL source;
  document provenance.
- Component licences: admixtools2 **GPL-3** (its DESCRIPTION); R, including
  `src/appl/lbfgsb.c` and `RNG.c`, **GPL-2-or-later**; L-BFGS-B 3.0 Fortran
  (Byrd, Lu, Nocedal, Zhu) **BSD-3**; MT19937 reference code (Matsumoto,
  Nishimura) **BSD-style**; the Broad C library shared by DATES, ADMIXTOOLS
  and EIGENSOFT (`DATES/src/nicksrc/LICENSE.txt`) **BSD-3**, compatible
  with MIT given attribution. To check: ALDER's, ROLLOFF's, READv2's,
  ADMIXTURE's and EIGENSOFT's own top-level licences.
- Not legal advice. A lawyer's review before 1.0.0 removes the residual risk.

Tasks:
- [x] **C1. Provenance audit.** *Done in 6fe4e9f.* For every module, record what it was written
  from (paper, documentation, black-box reference output, or reference
  source) and whether any GPL source was translated. Flag modules whose
  structure, names or comments follow admixtools2 R code or R's C source.
- [x] **C2. Rewrite flagged modules from neutral sources.** *Done: rewritten clean-room from docs/cleanroom/ (1a84440, 1c4cd69, f367488, a07192a, 4cea177, 77f3f64, 3c4984f, 734f710, 05663c2, f8412a5, 0a12912), or removed (589be68, 9ed980b); see PROVENANCE.md.* The audit's list (PROVENANCE.md): qpadm/fit.rs, lazadm.rs, qpgraph_fit.rs
  edge lengths, graph_search.py, qpfstats.py, graph_constraints.py,
  simulate.py, classic.py (admixtools2, GPL-3); all of alder/ (ALDER 1.03,
  academic use only); rcompat/lbfgsb.rs (R, GPL-2+; re-base on BSD 3.0);
  then the partial ones: graph_search.rs order, jackknife.rs, f4ratio.rs,
  resample.py, qpadm_models.py, rng.rs, dates/ld_decay.rs mesh helpers (DATES
  has no licence). Golden tests of rewritten modules should compare public
  outputs, not admixtools2's internal stages. L-BFGS-B from the
  paper and the BSD 3.0 Fortran (keeping R's 2.x behaviour, documented as
  behaviour); R's RNG from MT19937 (BSD) plus R's documented seeding,
  verified black-box against R; admixtools2 algorithms from Patterson et al.
  2012, Harney et al. 2021, Maier et al. 2023 and black-box outputs. The
  golden tests stay the oracle.
- [x] **C3. Remove source-line pointers into GPL code** *Done with the TASKS update after 6fe4e9f.* from comments and
  docs (file and line references), keeping behaviour descriptions.
- [x] **C4. NOTICE / THIRD-PARTY file** *Done in 6fe4e9f.* with the BSD attributions (Broad
  library, L-BFGS-B, MT19937) and any other required notices.
- [x] **C5. Process change (CLAUDE.md, `.cursor/rules`):** *Done in b556631.* implement from
  papers, documentation and reference outputs; do not translate or consult
  GPL source while writing code. The R source can still be run as an oracle.
- [x] **C6. Check the remaining tools' licences** *Done in 6fe4e9f.* (ALDER, ROLLOFF, READv2,
  ADMIXTURE, EIGENSOFT top-level) and record them in C4.
- [ ] **C7. Optional: legal review** of the result before tagging 1.0.0.

## D. Release blockers: dead code, docs, platform (P0)

- [x] **D1. Delete the unused "ergonomics" layer** *Done in 53c02fc.* (from b0ff796):
  `python/popmix/cli/{rich_output,formatters,statistics,batch}.py`,
  `python/popmix/reports/`, `python/popmix/multiqc/`, `python/popmix/schemas/`,
  their tests, and the `[cli]`, `[reports]`, `[multiqc]`, `[formats]` extras
  (or rebuild what is wanted properly with end-to-end tests). They crash on
  real results (`AttributeError: 'amplitude'`, `'r_squared'`), drop SEs via
  `hasattr(result, "std_errs")`, and `html_report.py:400-411` draws a
  made-up decay curve. AI [reviewer-verified]
  - includes the scientifically wrong helpers in `cli/statistics.py:381-415`
    (positive f3 "no admixture, high confidence"; Cohen's d; |Z| 2-3
    "significant"; own p-value approximation) PR, PD [reviewer-verified]
  - and the qpAdm "Good Fit"/PASS/FAIL labels that ignore feasibility, with
    the weight chart clipped to [0, 100]% (`reports/html_report.py:225-231,
    ~330`, `cli/formatters.py:165`, `cli/rich_output.py:372,421`) PR
- [x] **D2. README matches the program.** *Done in 619dc5a.* Remove or implement `qpadm-batch`,
  `dates-batch`, `--report`, `-o/--output-file`, yaml/parquet/ndjson formats,
  "Rich styling (default)", "schema-validated" JSON, `multiqc --module
  popmix` (README.md:654-760); `extract_f2(prefix)` without populations
  ("slow, computes all pairs") but `populations` is required; "streaming
  extraction for large datasets" true only for PACKEDANCESTRYMAP.
  AI, PD [confirmed for qpadm-batch, dates-batch, --report]
- [x] **D3. TESTING.md from real counts and commands** *Done in 240794c.* (claims 111/26/10 Rust
  tests and "26 passed"; documents `popmix f2 --snp-range`, which only
  `extract-f2` has; lists the dead modules). AI [reviewer-verified]
- [x] **D4. Test every documented command.** *Done in 619dc5a.* A test that runs each `popmix …`
  line in README/TESTING with `--help` and checks its flags exist. AI
- [x] **D5. Rename the 23 files with `:` in their names** *Done in 87ea173.*
  (`tests/golden/dates/Simulation_French_Yoruba/Simulation.out:1` …) so
  Windows can check out the repo; the release workflow's Windows job fails
  otherwise. AI [confirmed]
- [x] **D6. Remove or mark experimental the legacy `fit_graph`** *Done in 080e504.*
  (`qpgraph/fitting.rs`, Levenberg-Marquardt, two `#[allow(dead_code)]` at
  293 and 468, exported as `pm.fit_graph`, no golden test). AI, RS
- [x] **D7. Remove tracked scaffolding** *Done in fdec296.* (`.cursor/plans/ergonomics.md`),
  stray comments ("Actually, let's use…" at `f2_blocks.rs:446`). AI, RS
- [x] **D8. Fix the stale tolerance example** *Done in 0d0b1ec.* (hard-coded `0.0623`,
  `abs < 0.0001`) in CLAUDE.md and `.cursor/rules/004-testing.mdc`, which
  contradicts the tolerance table. AI

## E. Benchmarks against the reference tools, visible (P0)

- [x] **E1. Benchmark harness.** *Done in 05bb262.* For each analysis, run popmix and the
  reference tool on identical inputs, N repeats (default 5): wall time, CPU
  time, peak RSS, threads; median, IQR, coefficient of variation.
- [x] **E2. Agreement in the same run.** *Done in 05bb262.* Max absolute and relative difference
  of every output against the reference, pass/fail at the CLAUDE.md
  tolerances.
- [ ] **E3. Variation.** popmix across repeats and thread counts (1, 2, 4, 8,
  all), expected bit-identical; the reference tools' own run-to-run spread
  where they are stochastic (ADMIXTURE seeds, qpGraph starts, find_graphs).
  (popmix across 1/4/16 threads: done in 05bb262, identical in every case;
  the stochastic tools' spread remains, with their cases. ADMIXTURE: both
  tools' log-likelihood spread over seeds 1-3 is in its agreement.)
- [x] **E4. Coverage.** *Done: cases in `benchmarks/cases.py`.* admixtools2 (extract_f2, f2/f3/f4/qpdstat, qpadm,
  qpadm rotation, qpgraph, find_graphs, resampling), classic ADMIXTOOLS
  (qpDstat, qp3Pop, qpAdm, qpGraph), smartpca, ADMIXTURE, DATES, ROLLOFF,
  KING (plink2), READv2. qpfstats and ALDER are dropped: popmix 1.0 no
  longer has them (C2). Open findings from the trial runs:
  - classic qpAdm keeps another SNP set than admixtools2 and popmix
    (803,891 against 675,980 on the case's model), so weights differ by
    1e-3 and χ² by 1.5%: beyond qpAdm's printed precision, the case fails.
    Diagnosed: qpAdm keeps a SNP when every basis f4 is evaluable, an f4
    with a zero factor counting as 0 even where a population has no data,
    and keeps monomorphic SNPs; on a common SNP set its SEs and χ² still
    differ from admixtools2's. A difference of method, not reproduced
    (DEVIATIONS.md, "ADMIXTOOLS 7.0.2 (qpAdm)").
  - DATES with a pseudo-haploid source (French = Sardinian +
    Russia_Samara_EBA_Yamnaya, AADR v66.p1): popmix with
    `compat="dates-v753"` counts 0.06% more pairs per bin than DATES, and
    the fit moves from 152.9 to 153.8 generations. With diploid sources it
    matches; the case uses those until this is understood. Diagnosed:
    DATES counts every sample as diploid, and the Yamnaya source holds a
    diploid genome (I6729.DG) among pseudo-haploid samples; reproduced
    under the compat flag (DATES-0004), and the case uses Yamnaya again.
  - the stock `/home/a/DATES/bin/dates` aborts at random in mesh mode too
    (uninitialized `numbins`), and its `dowtjack` needs a missing
    libgsl.so.19: the case builds DATES with `dates_numbins_fix.patch` in
    the cache.
- [ ] **E5. Scaling curves** in SNPs, samples, populations and threads, up to
  full AADR v62 (17,629 samples × 1.23M SNPs).
- [ ] **E6. Report.** Results as JSON in the repo with tool versions and
  machine details; a script renders tables and charts to an HTML page
  (published as an artifact) and a summary in README. A small-data subset
  runs in CI so regressions show up.

## F. Performance and scalability (requirement)

Targets: at least as fast as each reference tool on every E benchmark;
full AADR analyses within a 62 GB machine; parallel speed-up close to
linear to 8+ threads where the work allows; results identical at every
thread count.

- [ ] **F1. Stream every format.** PLINK `.bed`, VCF and text EIGENSTRAT
  straight into allele frequencies, as packed `.geno` does; today
  `extract_f2` fallback (`lib.rs:1464`), `geno_fstats.rs:53`, `sfs.rs:78` call
  `read_all()` into `Array2<u8>` (21.7 GB for AADR). PD [reviewer-verified]
- [ ] **F2. DATES/ALDER/ROLLOFF read only the samples they need** from a
  prefix (`popmix dates` loads the whole matrix, `cli/__init__.py:1198`). PD
- [ ] **F3. Streaming `merge` and `convertf`** (row iterators; today both
  inputs plus output are in memory, ~65 GB for two AADR-sized files).
  PD, PG
- [ ] **F4. f2 extraction:** no copy of the frequency table (index through
  `keep`, `f2_blocks.rs:216-222`); one pass for f2 and FST; write block
  results in place; GEMM formulation (Σp² + Σq² − 2pq plus bias terms) for
  many populations; population-pair chunking with a memory limit
  (admixtools2's `maxmem`). Measured: 25 pops 13.4 s / 1.17 GB, 100 pops
  55.7 s / 3.48 GB; ~1,000 pops extrapolates to ~35 GB. PD
  [reviewer-verified]
- [ ] **F5. PCA eigensolver.** nalgebra `SymmetricEigen` is dense,
  single-threaded, O(n³), rerun every outlier iteration (n = 2000: 9.9 s;
  ~20 min per iteration at 10k). Use a parallel solver (faer) or a partial
  top-K Lanczos/randomized solver plus the trace for Tracy-Widom; bound the
  n×n accumulators per in-flight group (~800 MB each at n = 10k). PD
  [reviewer-verified]
- [ ] **F6. Release the GIL** in every long-running `lib.rs` function
  (`extract_f2`, `qpadm`, `qpadm_rotate`, `qpgraph`, `find_graphs`,
  `convertf`, `merge`, `read_*`); 0 `allow_threads` there today. RS
  [reviewer-verified]
- [ ] **F7. Thread pool that honours its contract.** The default of 4 applies
  only if `init_thread_pool` runs first; `get_threads()` creates the global
  pool (16 here) and `init_threads(8)` then returns `False`
  (`threading.rs:47`). Own a pool and `install` each entry point. RS
  [reviewer-verified]
- [ ] **F8. Fixed-order reductions everywhere,** with a thread-count
  determinism test for every parallel sum (A9 first).
- [ ] **F9. Graph search:** compute `precompute_f3` once per population set,
  not per candidate (`graph_search.rs:483-495`); `HashMap` population index
  instead of linear `pop_index` scans; port the Python evolutionary search
  (`python/popmix/graph_search.py`, 750 lines) to Rust. PD, RS, AI
- [ ] **F10. Zero-copy and fewer copies across the Python boundary:**
  `to_numpy` copies row by row (`lib.rs:188`); `#[pyo3(get)]` on `Vec` fields
  (129) clones on every access. RS
- [ ] **F11. Kernels:** check the genotype and LD kernels auto-vectorize
  (packed popcount, f2 inner loops); use safe SIMD (`pulp`/`wide`) or
  runtime-dispatched variants where it pays; avoid dense Kronecker
  temporaries (`qpadm/fit.rs:68,85`). RS, PD
- [ ] **F12. Benchmarks in-repo** with criterion (already a dev-dependency,
  unused) for the hot kernels, and profiling (flamegraphs) to choose work.
- [ ] **F13. Readers usable across threads:** map in `new()`, read via
  `&self` (`plink.rs:128`, `packedancestrymap.rs:143`), with `// SAFETY:`
  comments on the two `Mmap::map` calls. RS
- [ ] **F14. KING kinship speed:** popmix 3.9 s against PLINK 2's 0.84 s (16
  threads) for 220 samples x 1.2M SNPs in the E benchmark; packed popcount
  kernels as plink2 uses (F11).

## G. Scientific guardrails and honest reporting (P1)

- [ ] **G1. qpAdm verdict:** combine p, feasibility and SEs; no PASS/FAIL or
  "Good model fit" from p alone (`cli/__init__.py:974-975`); report the
  number of models tested in rotations and model sets (multiple testing,
  Harney et al. 2021). PR, PG
- [ ] **G2. qpWave wording:** rank is a lower bound limited by the right
  set; not rejecting rank 0 is not "Populations form a clade"
  (`qpwave.rs:103-108` hard-codes 0.05; `cli/__init__.py:1038-1040`). PR
- [ ] **G3. DATES/ALDER dates:** "years BP" is generations × 29 before the
  sample lived (`dates/dating.rs:254`, `cli/__init__.py:1254`): take the
  sample age, make generation time an option with uncertainty, and stop
  presenting R² as evidence of a single pulse. PR, AI
- [ ] **G4. Warn on few blocks or chromosomes:** leave-one-chromosome-out
  jackknife on 3 chromosomes; `f4_matrix.rs:120` accepts 2 blocks; qpAdm's
  χ² treats the jackknife covariance as known, inflated when (nl−1)(nr−1)
  nears the block count (document; consider a Hotelling-type correction).
  PR
- [ ] **G5. Warn when there is no genetic map** (2 Mb blocks,
  `fstats/blocks.rs:63`, shorter than 5 cM; SEs slightly small). PR
- [ ] **G6. Projected ancient samples:** default to lsqproject when
  projected samples have much missing data, or warn; implement shrinkage
  correction (PLAN.md A2 listed it; `pca/smartpca.rs:145`). PR
- [ ] **G7. Single pseudohaploid target in f3:** the CLI `f3` (f2 path) gives
  f3 = 0.19, z = 122, while `qp3pop -i` gives NaN with n = 0 (both as
  admixtools2); make them agree or warn. PG
- [ ] **G8. z = NaN when SE = 0,** not 0 (`f3_f4.rs:92-96`, f2 self-pairs).
  PD
- [ ] **G9. qpGraph and find_graphs:** report overfitting and
  identifiability (penalized score with an admixture-count penalty,
  out-of-sample score by default; `compare_fits` exists). PR
- [ ] **G10. DATES `min_source_count` of 1:** warn when sources have very
  few called samples. PR
- [ ] **G11. Kinship caveats:** KING −∞ for pseudohaploid samples, fixed READ
  thresholds; document coverage-dependent bias and contamination. PR
- [ ] **G12. Jackknife docstring:** the est/tot order in
  `fstats/jackknife.rs:186` is valid only for additive statistics (as
  admixtools2 `jack_pairarr_stats`); say so. PR [reviewer-verified algebra]
- [ ] **G13. Warn on the combinations users get wrong without noticing:**
  generalize the above into a warnings channel (Rust → Python → CLI), today
  almost absent.

## H. Practitioner features (P1)

- [ ] **H1. SNP filters** for extract-f2 and qpadm: transversions only
  (essential for non-UDG aDNA), a keep/exclude SNP list, minimum MAF. PG
- [ ] **H2. Sample selection on the CLI:** include/exclude individuals and
  relabel them (`inds`/`pops` exist only in Python); a population subset for
  `convertf`. PG
- [ ] **H3. Provenance in outputs:** f2 directory `metadata.json` and result
  files record all parameters (maxmiss, block size, auto_only, ploidy per
  sample, SNP count) and the popmix version. PG
- [ ] **H4. qpAdm CLI output complete:** SEs and z in text and TSV, rank
  tests, popdrop/nested models, SNP and block counts, rotation output in
  TSV; p-values and SEs in full precision or scientific notation (not
  `0.0000`, not `.6f`). PG, PR
- [ ] **H5. Strand-ambiguous SNPs kept on request** when both files share a
  strand (two AADR releases). PG
- [ ] **H6. qpadm `--compat`** on the CLI, as DEVIATIONS.md step 2 promises.
  AI
- [ ] **H7. Calibration by simulation:** a command giving the false-positive
  rate and power of a qpAdm or f4 design on msprime simulations. PR
- [ ] **H8. Admixture dating model checks:** single pulse vs extended
  admixture (residual test or two-exponential fit). PR
- [ ] **H9. Ascertainment awareness:** record the SNP panel in f2 metadata;
  warn when f2 directories from different panels are combined. PR
- [ ] **H10. Coverage-aware ancient-sample handling** across PCA, kinship and
  ploidy (G6, A6). PR

- [ ] **H11. qpAdm/qpWave on all SNPs of each f4** (admixtools2's
  `qpadm(prefix)`, which computes each f4 from the genotypes on the SNPs its
  populations share): an option beside today's f2-block path, which matches
  `qpadm(f2_from_geno(...))`. Found by the E benchmark (DEVIATIONS.md).

## I. Rust API and error design (P1)

- [ ] **I1. `PopmixError`:** `#[non_exhaustive]`; `PathBuf` paths with
  `#[source]` io errors (`io_with_path`, `error.rs:118`, drops the path except
  for NotFound; bare `?` in `npy.rs`); real variants instead of
  `Validation(String)` (31 uses) and `Numerical(String)` (40); remove the ten
  never-constructed variants; NPY format errors as `InvalidFileFormat`. RS
  [reviewer-verified]
- [ ] **I2. Typed Python exceptions:** `FileNotFoundError`, `OSError`,
  `KeyError`/`ValueError` for populations and validation, popmix exception
  subclasses; `From<PopmixError> for PyErr` so `.map_err(to_py_err)` becomes
  `?` (everything is `RuntimeError` today, `lib.rs:36`). RS, PD
- [ ] **I3. Deliberate public API:** private fields with accessors on
  `F2Blocks`/`AlleleFrequencies`; `#[non_exhaustive]` options and results;
  `&Options` not `Option<Options>`; consistent `impl AsRef<Path>`; `impl
  FromStr for GenoFormat`; `rcompat` and the second graph fitter
  `pub(crate)` or documented; distinct names for dates and qpgraph
  `FitOptions`; decide whether ndarray/nalgebra types are public. RS
- [ ] **I4. Move the extraction pipeline out of the bindings:** format
  detection, streaming choice, ploidy detection and relabelling into
  `popmix_io::extract_f2(prefix, &ExtractOptions)` behind a `GenotypeSource`
  trait; `lib.rs` (3,153 lines) becomes thin wrappers, split by domain;
  delegate `best_clade`/`worst_pair` (`lib.rs:2079-2083,2345-2349`) to core.
  RS, AI
- [ ] **I5. One linear-algebra library** in public types and hot paths. RS
- [ ] **I6. Don't promise a stable Rust API at 1.0** until I1-I4 are done
  (Python package 1.0 is separate). RS

## J. Python API consistency (P2)

- [ ] **J1. One input convention:** `dates()` takes a matrix while `alder`,
  `smartpca`, `king` take a prefix; `qpadm` `sources` vs `qpwave` `left`. PD
- [ ] **J2. One result convention:** dicts (`f2`), pyclasses, DataFrames with
  admixtools column names (`est`/`se`); counts as ints, not floats
  (`n_snps: 21077.0`). PD, AI
- [ ] **J3. Type hints and NumPy-format docstrings** on the raw PyO3
  functions (CLAUDE.md rule). PD
- [ ] **J4. pyclass `module="popmix._core"`, names, `frozen` where possible,
  pickling** (reprs show `builtins.PyQpAdmRotateResult`). RS, AI
- [ ] **J5. `f4` help text** says "D-statistic / ABBA-BABA" but prints f4. PG

## K. Tests, CI and process (P1)

- [ ] **K1. CI test workflow:** `cargo test`, `pytest tests/python`, and the
  golden tests that ship their data (family_chr20_22, rlang, family_4pop) on
  Linux, macOS and Windows, x86_64 and aarch64; a skipped golden fails CI
  (`POPMIX_REQUIRE_GOLDEN=1`). AI, PD
- [ ] **K2. Portable AADR goldens:** a Zenodo-hosted AADR subset with
  checksums; nightly or self-hosted run of the AADR goldens; record
  `sessionInfo()` and tool versions in every reference. PD
- [ ] **K3. Cross-format reference tests:** the same individuals through VCF,
  PLINK, PED, text and packed EIGENSTRAT give identical f2, including a
  pseudohaploid population; VCF → EIGENSTRAT/PLINK round-trips checked
  against plink2 and EIGENSOFT convertf; merge checked against mergeit and
  `plink --bmerge`. PG (Merge against mergeit done with A10: plink 1.9
  --bmerge stops at the first SNP needing a strand flip, so it serves only
  as a list of SNPs to `--flip`.)
- [ ] **K4. Repeated-population and degenerate-model goldens:**
  f4(A,B;A,C) = f3(A;B,C), f3(A;A,B) = 0 as properties; dof 0; overlapping
  left/right; against admixtools2. PD, PR
- [ ] **K5. Replace mock-based and can't-fail tests:** `MockF2Result`,
  `MockQpAdmResult` (`test_statistical_rigor.py:13-35`), single-character
  and `or`-chained assertions (`test_rich_output.py:276-279,422`,
  `test_cli.py:61`), the mismatched docstring at `test_statistical_rigor.py:101`,
  the tag-presence checks in `test_html_reports.py`, the weights-sum-to-1
  qpadm test (`qpadm_impl.rs:227-269`). Most go with D1. AI, PD
- [ ] **K6. Mutation testing** (`cargo-mutants`, `mutmut`) on core statistics
  and Python helpers, to find tests that cannot fail. AI
- [ ] **K7. Process rules in CLAUDE.md:** a feature is done only when it is
  reachable from the CLI or `popmix.__all__` and an end-to-end test runs it
  on a real result object; no mock result classes when the real type can be
  built; no `or`-chained assertions. AI
- [ ] **K8. Tests that pin the release workflow by running it:** a test tag
  run of `release.yml` on every platform before 1.0.0. AI

## L. Packaging and release (after A-E)

- [ ] **L1. popmix-io/core publishable:** `version` on the path dependency,
  `publish = false` on popmix-python, `rust-version` on every crate. RS
- [ ] **L2. Unused dependencies:** `ordered-float`; `proptest` and
  `criterion` unless F12/B6 use them. RS
- [ ] **L3. PyO3 upgrade** (0.22 `*_bound` API deprecated from 0.23); abi3
  if rust-numpy allows. RS
- [ ] **L4. Release:** create `trianglegrrl/popmix`, push, run `release.yml`
  by hand, tag `v1.0.0` (PyPI trusted publisher is configured).
- [ ] **L5. Bioconda:** recipe from the PyPI sdist and sha256, local build,
  PR to bioconda-recipes.
- [ ] **L6. nf-core:** pin all 21 modules to the released 1.0.0 container
  (they pin `popmix=0.1.0`, `biocontainers/popmix:0.1.0--py39h2e5f10c_0`);
  add a merge/convertf test from VCF; test data to nf-core/test-datasets;
  modules PR. PG

## M. Minor cleanups (P2)

- [ ] **M1. Duplicated logic:** three confidence-interval implementations
  (`formatters.py:34`, `rich_output.py:27`, `statistics.py:26`; go with D1).
  AI
- [ ] **M2. Idiom:** `f64::midpoint` (9 places), `let … else` instead of
  `unwrap` after `is_none` (`ped.rs:184`, `graph.rs:240`), worthwhile
  `clippy::pedantic` lints (1,782 warnings, mostly noise). RS
- [ ] **M3. Document why graph search is split across Rust and Python** (or
  resolve it with F9). AI

## N. Agent skill for popmix (P1)

- [ ] **N1. Skill layout with progressive disclosure:** `skills/popmix/SKILL.md`
  (what popmix does, when to use which command, input formats, an index of
  commands) plus one file per CLI command, `skills/popmix/commands/<command>.md`
  (subcommands get their own files), each complete: synopsis, every
  argument and option with type, default and meaning, input and output
  files, worked examples on real data, how to read the result, caveats, and
  the reference tool and function it corresponds to.
- [ ] **N2. 1:1 mapping enforced by a test:** walk the Click command tree;
  every command has exactly one doc file and every doc file one command;
  every option and argument of the command appears in its file with the
  right default; fails when the CLI and the skill drift.
- [ ] **N3. Examples that run:** the examples in each command file execute in
  the test suite (small bundled data), so they cannot go stale.
- [ ] **N4. Ship it:** include the skill in the package (installed with
  popmix) with a command or documented step to install it for Claude
  Code (`~/.claude/skills/popmix`); version it with popmix.

## O. AADR v66.p1 as the reference dataset (P0)

- [x] **O1. Read TGENO** *Done in 66b3375.* (AdmixTools `transpose_packed`, one row per
  individual) in popmix, streaming; it suits reading a subset of
  individuals. v66.p1 1240K `.geno` begins `TGENO 23089 1233013`.
- [x] **O2. TGENO misdetected as text EIGENSTRAT.** *Done in 66b3375.* `dataset_info` on the v66.p1
  prefix reports `format: eigenstrat`; analyses would read garbage. Detect
  the header and support it (O1) or refuse with a clear error. [confirmed]
- [x] **O3. Build AdmixTools ≥ 8 convertf** *Done in 26e9dbb.* as a reference tool (the system
  EIGENSOFT convertf fails: "bad ID (SNP): TGENO ..."); record its version.
- [ ] **O4. Migrate every AADR golden from v62.0 to v66.p1** (`/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB`):
  v66.p1 population labels carry no data-type suffix (`French`, not
  `French.DG`; the suffix is on sample IDs such as `Loschbour.AG`), so
  each reference names its individuals or relabels by suffix. Regenerate
  all references and record the dataset checksums
  (`CHECKSUMS.normalized.md5`).
  *References regenerated (tests/golden/aadr/build_aadr_v66p1.py builds the
  labelled dataset). Open: three live checks fail on v66.p1 data and are
  not papered over: DATES `aadr_french_q10-default_result` (DATES-0003
  moves the jackknife by 0.16, beyond the check's bound),
  `ancient_5pop-qpgraph_fit` (L-BFGS-B ends 1.9e-6 above admixtools2), and
  two 2-admixture graphs at admixtools2's weights (1.6e-10 apart).*
- [x] **O5. Update CLAUDE.md, README, TESTING.md and generator defaults** to
  v66.p1 paths; keep v62.0 only where a test pins a specific old result.
- [x] **O6. Data ethics check in CI:** *Done in 803ae88.* a test that fails if any individual or
  population removed in a `.p1` release appears in the repository or test
  data.

## P. `popmix convertf`, a full convertf replacement (P1)

- [ ] **P1. Every EIGENSOFT/AdmixTools format in and out:** PACKEDANCESTRYMAP,
  TGENO (transpose_packed), EIGENSTRAT, ANCESTRYMAP, PACKEDPED, PED, plus
  VCF in (and VCF out). Today TGENO is unsupported in both directions and
  VCF input is inverted (A1). (VCF input fixed in A1 and A7, TGENO in and
  out in 66b3375; VCF out remains.)
- [ ] **P2. Parameter-file compatibility:** `popmix convertf -p file.par` with
  convertf's keys (genotypename, snpname, indivname, outputformat,
  genooutfilename, snpoutfilename, indoutfilename, poplistname, ...), so
  pipelines can swap it in.
- [ ] **P3. Subsetting and relabelling:** population list, individual
  include/exclude, SNP list, chromosome and region, relabel populations
  (with H2).
- [ ] **P4. Fast:** streaming, parallel, bounded memory; benchmarked against
  convertf in section E.
- [ ] **P5. Checked against convertf on every format pair** (golden), with
  K3.
- [x] **P6. Packed headers carry AdmixTools' hashes.** *Done in 281226d.* popmix wrote
  `GENO <n> <m> 0`, so AdmixTools 8 convertf aborts on popmix's packed files
  ("indiv file has changed since genotype file was created"). Found by
  K3's test. [confirmed]

