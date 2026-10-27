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
- The 1.0.0 release waits for **every task in this file** (2026-10-27),
  C7 (legal review) alone excepted; L4-L6 are the release itself, done
  last. Order: "Path to 1.0" below.
- Reference data is **AADR v66.p1** from now on (section O). It removes 161
  modern individuals of Jacobs et al. 2019 (Cell 177:1010) that v62.0 and
  v66.0 contained; none of them, nor their 27 populations, appear in the
  repository, the reference cache or the nf-core test data (checked
  2026-10-25).
- popmix ships an agent skill documenting every CLI command (section N).
- `popmix convertf` becomes a full replacement for EIGENSOFT/AdmixTools
  convertf (section P).
- 1.0 is released for **Linux only** (2026-10-27): K1's macOS and Windows
  jobs run once the repository is public, after 1.0.
- H7 takes simulated genotypes as input (users run msprime themselves;
  popmix never imports it, msprime being GPL-3) and reports false-positive
  rate and power; the README's caveats section says so.
- K2's AADR subset goes on Zenodo (the user's account); K8's dry run uses a
  TestPyPI trusted publisher the user sets up. C7 is with the user.
- **Small source files** (2026-10-27): split by concern, no source file past
  about 800 lines, long test modules in their own files (I7).

## Path to 1.0

Work in this order; within a phase, top to bottom. Each phase makes the next
cheaper: the refactors come first because every later task edits the code
they move, and the public API is settled before features add to it.

1. **Foundations.** K7 (process rules), L3 (PyO3 upgrade), L2, L1, I1, I2,
   I4 (pipeline out of the bindings, `lib.rs` split), I7 (small files), I3,
   I5, M2, M1.
2. **The API 1.0 freezes.** G13 (warnings channel), J1, J2, J4, J3, J5, M4,
   M5, H6, I6, M3.
3. **Honest science.** G8, G12, G2, G1, G3, G4, G5, G7, G10, G11, G6, G14,
   G9.
4. **Practitioner features.** H3, H9, H1, H2 with P3, H4, H5, H11, P1 (VCF
   out), P2, H8, H10, H7.
5. **Performance left over.** F4 (pair chunks past ~1,000 populations, on
   I4's pipeline), F14 (KING at 3,970 samples).
6. **Tests and CI.** K5, K4, K3 with P5, P4, K6, K1 (Linux; macOS and
   Windows after 1.0), K2, K8.
7. **Release.** L4, L5, L6. C7 stays open.

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
- [x] **A11. qpWave rank when every test is rejected:** `qpwave.rs:113`
  falls back to the highest rank tested (`unwrap_or(max_rank)`), so with
  every rank rejected it reports one rank too few: `--left French_1,French_2`
  says rank 0 though rank 0 has p = 0.039. All rejected means rank above
  `max_rank`; report `max_rank + 1`. Check qpAdm's feasibility uses the same.
- [x] **A12. `extract-f2 --snp-range` without `--populations`** finds 0
  populations, writes an empty f2 directory and exits 0; with
  `--snp-range`, `--maxmiss`, `--auto-only` and `--pseudohaploid` are
  ignored. Error out, or honour them.

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
- [x] **B14. CLI options ignored without a word:** `qp3pop --outgroup-mode`
  and `--no-bias-correction` with `--f2-dir`; `sfs --joint` with
  `--projection`, `--folded`, `--all-chromosomes`. Error out when they
  cannot apply.
- [x] **B15. Broken pipe:** `popmix kinship ... | head` prints
  `Error: [Errno 32] Broken pipe` and exits 120; exit quietly.
- [x] **B16. Sex lost on conversion:** `convertf` to PLINK writes sex 0 and
  `merge` writes `U` for every sample although the `.ind` has M/F.

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
- [x] **E3. Variation.** popmix across repeats and thread counts (1, 2, 4, 8,
  all), expected bit-identical; the reference tools' own run-to-run spread
  where they are stochastic (ADMIXTURE seeds, qpGraph starts, find_graphs).
  (popmix across 1/4/16 threads: done in 05bb262, identical in every case;
  ADMIXTURE, qpGraph and find_graphs: each tool's spread over seeds 1-3 is
  in its case's agreement; qpGraph's best of 10 starts spreads by 1e-13
  (popmix) and 0 (admixtools2), find_graphs' evolutionary search by 8e-12.)
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
- [x] **E5. Scaling curves** in SNPs, samples, populations and threads, up to
  full AADR v66.p1 (23,089 samples × 1.23M SNPs). *Done: `benchmarks/scaling.py`,
  results in 5c1c90c; the findings are F14, F17 and F18.*
- [x] **E6. Report.** Results as JSON in the repo with tool versions and
  machine details; a script renders tables and charts to an HTML page
  (published as an artifact) and a summary in README. A small-data subset
  runs in CI so regressions show up. *Done: report_html.py, README
  (1919738), ci.yml with benchmarks/smoke.py; the workflow has not yet run
  on GitHub.*

## F. Performance and scalability (requirement)

Targets: at least as fast as each reference tool on every E benchmark;
full AADR analyses within a 62 GB machine; parallel speed-up close to
linear to 8+ threads where the work allows; results identical at every
thread count.

- [x] **F1. Stream every format.** PLINK `.bed`, VCF and text EIGENSTRAT
  straight into allele frequencies, as packed `.geno` does; today
  `extract_f2` fallback (`lib.rs:1464`), `geno_fstats.rs:53`, `sfs.rs:78` call
  `read_all()` into `Array2<u8>` (21.7 GB for AADR). PD [reviewer-verified]
  Progress: PLINK .bed is read in place like packed .geno (memory-mapped SNP-major reader: rows, ploidy, parallel frequencies) by extract_f2, the SFS, KING/READ, PCA, ADMIXTURE and the genotype f-statistics (5ef0734): KING at 1,000 samples 17.7 s / 2.1 GB to 13.3 s / 1.0 GB. Text EIGENSTRAT whose lines are one character per sample (as popmix, convertf and admixtools2 write it) is read in place too, by every analysis above and by convertf and merge: extract_f2 of 16 populations from 473 x 1.2M text 6.9 s / 2.0 GB to 0.9 s / 1.2 GB; other text (CRLF, padded) goes through the reader. VCF is read once (it was read twice), lines parsed in parallel batches while a thread decompresses, genotypes held packed (2 bits, not a byte) and read through the same SnpMajor code: 1000G chr22 (2,504 x 1.06M) to packed 283 s / 2.8 GB to 16-22 s / 0.8 GB, identical output. Done: every format is read in place or, for VCF, packed in memory.
- [x] **F2. DATES/ALDER/ROLLOFF read only the samples they need** from a
  prefix (`popmix dates` loads the whole matrix, `cli/__init__.py:1198`). PD
  Done: `read_genotypes(prefix, populations=[...])` reads only those samples; `popmix dates` uses it; ROLLOFF already streamed (9f0a1ac).
- [x] **F3. Streaming `merge` and `convertf`** (row iterators; today both
  inputs plus output are in memory, ~65 GB for two AADR-sized files).
  PD, PG
  Progress: convertf from packed .geno or .bed to EIGENSTRAT, ANCESTRYMAP, PACKEDANCESTRYMAP or PLINK streams one SNP at a time (0bcf74c): 3,970 samples 15 s, 1.3 GB resident (mostly the mapped input) instead of the 4.9 GB matrix. `merge` of two packed .geno or .bed inputs streams too (MergedRows, read in place): 2 x 1,201 samples 22.7 s / 6.3 GB to 8.7 s / 1.3 GB, byte-identical. PED and TGENO output (stored sample by sample) are written from 256 MB blocks of samples read in place, from convertf and merge (SampleGenotypes): AADR 1,201 samples to TGENO 17.2 s / 2.1 GB to 2.6 s / 1.05 GB, to PED 79 s / 2.1 GB to 18 s / 1.0 GB (5.9 GB written), byte-identical. Done.
- [ ] **F4. f2 extraction:** no copy of the frequency table (index through
  `keep`, `f2_blocks.rs:216-222`); one pass for f2 and FST; write block
  results in place; GEMM formulation (Σp² + Σq² − 2pq plus bias terms) for
  many populations; population-pair chunking with a memory limit
  (admixtools2's `maxmem`). Measured: 25 pops 13.4 s / 1.17 GB, 100 pops
  55.7 s / 3.48 GB; ~1,000 pops extrapolates to ~35 GB. PD
  [reviewer-verified]
  Progress: block sums are matrix products and the kept SNPs are indexed in place (10dd2f5): `--maxmiss 1`, 4 threads, 128 populations 47.7 s to 15.7 s, 497 populations 676 s / 19.8 GB to 158 s / 13.8 GB; no block value differs from admixtools2 (16-256 populations). Then a chromosome at a time (F2Blocks::from_parts; blocks never span chromosomes, the map-or-base-pairs choice made from every SNP and checked against the kept ones at the end), block sums written in place: 497 populations 8.5 GB to 2.2 GB (1.5 GB of it the mapped input), with --maxmiss 1 13.8 GB / 158 s to 8.8 GB / 121 s, the same bytes. Left: past about 1,000 populations the blocks themselves (pairs x blocks) dominate; writing them to disk by pair chunks, as admixtools2's maxmem, would bound that.
- [x] **F5. PCA eigensolver.** nalgebra `SymmetricEigen` is dense,
  single-threaded, O(n³), rerun every outlier iteration (n = 2000: 9.9 s;
  ~20 min per iteration at 10k). Use a parallel solver (faer) or a partial
  top-K Lanczos/randomized solver plus the trace for Tracy-Widom; bound the
  n×n accumulators per in-flight group (~800 MB each at n = 10k). PD
  [reviewer-verified]
  Progress: faer's sequential eigensolver (5x nalgebra's at 2,000 samples), 27f0ab6. A partial top-K solver is not needed at the sizes measured. Done: at 2,000 samples the eigensolver takes 2.8 s of 115 s; X X' is accumulated in disjoint tiles of its upper triangle (n^2/2 values in all), not one n x n matrix per thread.
- [x] **F6. Release the GIL** in every long-running `lib.rs` function
  (`extract_f2`, `qpadm`, `qpadm_rotate`, `qpgraph`, `find_graphs`,
  `convertf`, `merge`, `read_*`); 0 `allow_threads` there today. RS
  [reviewer-verified]
  Done: every long-running binding runs through `compute` (GIL released, popmix's pool installed), d7be64f.
- [x] **F7. Thread pool that honours its contract.** The default of 4 applies
  only if `init_thread_pool` runs first; `get_threads()` creates the global
  pool (16 here) and `init_threads(8)` then returns `False`
  (`threading.rs:47`). Own a pool and `install` each entry point. RS
  [reviewer-verified]
  Done: popmix owns its pool (4 threads until `init_threads(n)`, which works at any time), d7be64f.
- [x] **F8. Fixed-order reductions everywhere,** with a thread-count
  determinism test for every parallel sum (A9 first).
  Progress: thread-count determinism tests for PCA, KING, f2 blocks, the frequency pass and the genotype f-statistics (all byte-identical in the benchmarks). Done: tests/python/test_thread_determinism.py runs PCA, KING, READ, ADMIXTURE, ROLLOFF, DATES, f2 blocks, convertf to every streamed format and merge at 1, 4 and 7 threads and requires the same bytes.
- [x] **F9. Graph search:** compute `precompute_f3` once per population set,
  not per candidate (`graph_search.rs:483-495`); `HashMap` population index
  instead of linear `pop_index` scans; port the Python evolutionary search
  (`python/popmix/graph_search.py`, 750 lines) to Rust. PD, RS, AI
  Progress: measured: caching precompute_f3 across the 1,965 graphs of the find_graphs benchmark changed nothing (1.3-1.5 s with or without), so it was not kept; the fits' optimization is the cost. The evolutionary search is out of 1.0 (C2). Done: popmix's exhaustive search runs the benchmark in 2.4 s against admixtools2's 45 s.
- [x] **F10. Zero-copy and fewer copies across the Python boundary:**
  `to_numpy` copies row by row (`lib.rs:188`); `#[pyo3(get)]` on `Vec` fields
  (129) clones on every access. RS
  Progress: `GenotypeMatrix.to_numpy` is one copy (f14f44e; 1.2M x 473 in 0.13 s). The Vec getters are on small result objects that the Python wrappers read once, so their clones were left. Done.
- [x] **F11. Kernels:** check the genotype and LD kernels auto-vectorize
  (packed popcount, f2 inner loops); use safe SIMD (`pulp`/`wide`) or
  runtime-dispatched variants where it pays; avoid dense Kronecker
  temporaries (`qpadm/fit.rs:68,85`). RS, PD
  Progress: KING's pair kernel has an AVX2/POPCNT build chosen at run time (2.6x; wheels target baseline x86-64); f2 blocks, PCA and the projections are matrix products; ADMIXTURE skips noise-level steps; .snp/.bim parse in parallel and the genotype f-statistics sum blocks in parallel (0454e55). Done: the KING kernel compiles to AVX2 popcount (vpshufb); ADMIXTURE's loops are specialized on K; qpAdm's Kronecker design matrix is at most about 150 x 30, milliseconds, and was left.
- [x] **F12. Benchmarks in-repo** with criterion (already a dev-dependency,
  unused) for the hot kernels, and profiling (flamegraphs) to choose work.
  Done: `crates/popmix-core/benches/kernels.rs` (criterion): PCA, KING, f2 blocks, ADMIXTURE on synthetic data; CI builds them (5531289).
- [x] **F13. Readers usable across threads:** map in `new()`, read via
  `&self` (`plink.rs:128`, `packedancestrymap.rs:143`), with `// SAFETY:`
  comments on the two `Mmap::map` calls. RS
  Done: the PLINK and packed readers map once (OnceLock) and read through `&self`, with SAFETY comments (9194e5e).
- [ ] **F14. KING kinship speed:** popmix 3.9 s against PLINK 2's 1.2 s (16
  threads) for 220 samples x 1.2M SNPs in the E benchmark; packed popcount
  kernels as plink2 uses (F11). At 3,970 samples 547 s and 7.4 GB against
  PLINK 2's 201 s. Threads, as in PCA (F17): at 2,000 samples 215 s on 1
  thread, 131 s on 4, 158 s on 16, CPU time 216 s to 1,527 s.
  Progress: tiled pair counts, AVX2/POPCNT kernel chosen at run time, words built per sample (911385b): 1,000 samples 61 s to 36 s (1 thread), 38 s to 18 s (4); pair counting 45 s to 20 s and 24 s to 7.6 s. Rerun (2026-10-26): E benchmark 3.9 s to 2.1 s (PLINK 2 1.2 s); 3,970 samples 547 s to 216 s and 3.8 GB (PLINK 2 201 s); 2,000 samples 105 s on 1 thread, 52 s on 4 (PLINK 2 135 s at its default). Then the E case's fixed costs: rows decoded in parallel and words built for all samples in one pass (read 250 to 190 ms), the SNP list moved rather than copied (open 470 to 160 ms), and the KING command writes its table without pandas (0.3 s of import): 1.5-1.9 s to 0.95 s at 16 threads, PLINK 2 0.78 s, same output. Then .bim and .snp files parsed in parallel pieces (parse_lines; 160-240 ms to well under 100 ms for 1.2M lines, which every command pays): KING 0.78-0.96 s, level with PLINK 2's 0.78 s; extract_f2 of 16 populations 1.01 s to 0.68 s. The pair kernel is AVX2 popcount (vpshufb) already; blocking four pairs per load gained nothing. At 3,970 samples popmix is within 7% of PLINK 2 (216 s against 201 s).
- [x] **F15. CLI startup:** `import popmix` loads pandas eagerly (through
  `popmix.admixture`), about a second on every command; import it lazily.
  Done: submodules load on first use; `popmix --help` 0.33 s to 0.04 s, 8d5ae2b.
- [x] **F16. ADMIXTURE speed:** popmix 19.3 s against ADMIXTURE's 13.9 s in
  the E benchmark.
  Progress: no passes for settled samples and SNPs, P-step log-likelihood reused (577c8fc): 53.0 s to 21.9 s (1 thread), 22.7 s to 6.3 s (4); ADMIXTURE 13.6 s on one thread. Rerun (2026-10-26): E benchmark 6.8 s on 4 threads (ADMIXTURE 13.9 s), 16.7 s on 1. Then one logarithm per genotype (ln a^h b^(2-h), chosen without a branch), reciprocals for the derivative weights, per-SNP and per-sample products shared, the inner loops specialized on K (2-8): one thread 15.7 s to 10.0 s against ADMIXTURE's 14.2 s, the same log-likelihood, identical at 1, 4 and 16 threads. (One-thread runs vary with the core they land on: 84 ms for a Q-step's sums on a P core, 150 ms on an E core.) Done.
- [x] **F17. PCA cost grows with the square of the samples:** in the sample
  scaling curve popmix takes about 4x per doubling (5.8 s at 250 samples,
  62.7 s at 1,000, 225 s at 2,000, 1,079 s at 3,970) and smartpca about 2x
  (63, 174, 406, 1,385 s); at 2,000 samples popmix uses 2.3x smartpca's CPU
  time. Threads do not help either: at 2,000 samples popmix takes 362 s on 1
  thread, 223 s on 4 and 286 s on 16, while its CPU time grows from 363 s to
  3,029 s, so the extra threads spin or contend instead of working. Profile
  the PCA path (the samples x samples product, its eigensolver, the parallel
  split), and make it scale like smartpca's with samples and threads.
  Progress: tiled upper-triangle X X', panel products for loadings and projections (27f0ab6): 1,000 samples 100 s to 65 s (1 thread), 53 s to 42 s (4). The GEMM runs at about 30 GFLOP/s per core. On this machine (i7-1360P, 4 P + 8 E cores, 16 logical CPUs, about 1.5 of them busy with other processes) FLOP-bound work stops scaling at 4-8 threads, and at 14-16 threads rayon workers are descheduled while the others spin, so 16 threads is slower than 4 for every kernel: the CPU-time growth in the scaling curves comes from that, not from popmix. Rerun (2026-10-26): 5.9, 14.3, 38.3, 143, 516 s at 250 to 3,970 samples (smartpca 63, 94, 174, 406, 1,385 s); at 2,000 samples 213 s on 1 thread, 119 s on 4, 136 s on 16. Faster than smartpca at every size, but its lead shrinks with samples (10.7x at 250, 2.7x at 3,970): the n^2 product still dominates. Measured (2026-10-27, 2,000 samples, 4 threads): X X' 97-108 s of 115 s (decoding 3 s, eigensolver 2.8 s), 49 GFLOP/s. That is this laptop's sustained ceiling: the same tiles run 85 GFLOP/s on 4 P cores for the first seconds and settle at 45 when its power limit sets in (OpenBLAS: 42 GFLOP/s on one P core, 18 on an E core). faer's GEMM in place of matrixmultiply gave nothing sustained. Both tools cost n^2 m; smartpca's lead-shrinking curve is its fixed costs dominating at few samples. Done: what is left is hardware, or a different method (randomized top-K), which would no longer give smartpca's full spectrum.
- [x] **F18. `extract_f2` runs on one core:** 128 populations take 16 s at
  every thread count from 1 to 16, CPU time equal to wall time. Parallelize
  the per-SNP frequencies and the per-block pair sums (populations x
  populations x blocks), which dominate past 64 populations: 51.6 s and
  6.9 GB at 497.
  Done: frequencies of every SNP in parallel from a mapped file, integer counts: 128 populations 16.4 s to 7.4 s (1 thread) and 3.9 s (4), byte-identical, ff63790.

## G. Scientific guardrails and honest reporting (P1)

- [x] **G1. qpAdm verdict:** combine p, feasibility and SEs; no PASS/FAIL or
  "Good model fit" from p alone (`cli/__init__.py:974-975`); report the
  number of models tested in rotations and model sets (multiple testing,
  Harney et al. 2021). PR, PG
  *Done in 4226c29:* `QpAdmResult.assess(alpha)` (status not tested /
  rejected / infeasible / not rejected, notes for weights within 2 SE of 0
  and simpler popdrop models that also fit); the CLI prints it, takes
  `--alpha`, JSON `assessment`; rotations and model sets warn with the
  number of models and set `attrs["n_models"]`.
- [x] **G2. qpWave wording:** rank is a lower bound limited by the right
  set; not rejecting rank 0 is not "Populations form a clade"
  (`qpwave.rs:103-108` hard-codes 0.05; `cli/__init__.py:1038-1040`). PR
  *Done in 43da0b9:* `alpha` option (Rust, Python, `--alpha`), in the result
  and JSON; `is_clade` removed; the CLI calls the rank a lower bound and
  says rank 0 not rejected is not a clade.
- [x] **G3. DATES/ALDER dates:** "years BP" is generations × 29 before the
  sample lived (`dates/dating.rs:254`, `cli/__init__.py:1254`): take the
  sample age, make generation time an option with uncertainty, and stop
  presenting R² as evidence of a single pulse. PR, AI
  *Done in 0e72ab3:* `sample_age(_se)` and `generation_time_se` (Rust,
  Python, CLI) enter the years and their SE; the CLI says "before the
  admixed samples lived" without a sample age, and that R² is not evidence
  of a single pulse.
- [x] **G4. Warn on few blocks or chromosomes:** leave-one-chromosome-out
  jackknife on 3 chromosomes; `f4_matrix.rs:120` accepts 2 blocks; qpAdm's
  χ² treats the jackknife covariance as known, inflated when (nl−1)(nr−1)
  nears the block count (document; consider a Hotelling-type correction).
  PR
  *Done in 4d7114c:* warnings below 20 f2 blocks, for qpAdm/qpWave
  chi-squared inflated past 10% by (g − 1)/(g − p − 2) (Hotelling; or
  g <= p + 2), and for DATES/ROLLOFF jackknives over fewer than 10
  chromosomes; documented, no correction applied.
- [x] **G5. Warn when there is no genetic map** (2 Mb blocks,
  `fstats/blocks.rs:63`, shorter than 5 cM; SEs slightly small). PR
  *Done in afdd936:* `Warning::NoGeneticMap` from extract_f2, genotype-path
  statistics and f2_from_frequencies.
- [ ] **G6. Projected ancient samples:** default to lsqproject when
  projected samples have much missing data, or warn; implement shrinkage
  correction (PLAN.md A2 listed it; `pca/smartpca.rs:145`). PR
- [ ] **G7. Single pseudohaploid target in f3:** the CLI `f3` (f2 path) gives
  f3 = 0.19, z = 122, while `qp3pop -i` gives NaN with n = 0 (both as
  admixtools2); make them agree or warn. PG
- [x] **G8. z = NaN when SE = 0,** not 0 (`f3_f4.rs:92-96`, f2 self-pairs).
  PD
  *Done in bbc1e33:* f3/f4, genotype-path qpdstat/qp3pop, f4ratio and
  `compare_fits` give z and p NaN when the SE is 0; f4 with pop1 = pop2
  or pop3 = pop4 is exactly 0 (admixtools2: rounding noise and a z from
  it; ADMIXTOOLS2-0008).
- [ ] **G9. qpGraph and find_graphs:** report overfitting and
  identifiability (penalized score with an admixture-count penalty,
  out-of-sample score by default; `compare_fits` exists). PR
- [ ] **G10. DATES `min_source_count` of 1:** warn when sources have very
  few called samples. PR
- [ ] **G11. Kinship caveats:** KING −∞ for pseudohaploid samples, fixed READ
  thresholds; document coverage-dependent bias and contamination. PR
- [x] **G12. Jackknife docstring:** the est/tot order in
  `fstats/jackknife.rs:186` is valid only for additive statistics (as
  admixtools2 `jack_pairarr_stats`); say so. PR [reviewer-verified algebra]
  *Done in c2d8285:* `jackknife_from_loo` and `jackknife_cov_from_loo` say
  their rebuilt total is the full-data value only for additive statistics
  (for D it is admixtools2's mean of leave-one-out ratios); a test pins it.
- [x] **G13. Warn on the combinations users get wrong without noticing:**
  generalize the above into a warnings channel (Rust → Python → CLI), today
  almost absent.
  *Done in 5baf897 (the channel; the warnings themselves are G4-G11):*
  `popmix_core::Warning` and a `warn: &mut dyn FnMut(Warning)` sink in Rust
  (since e309932); in Python each is a `popmix.PopmixWarning` (a
  UserWarning, filterable by class); the CLI prints each as one
  `Warning: ...` line on stderr. The CLI's own population-check notice
  goes the same way.
- [ ] **G14. qpGraph past a zero-bound kink:** where a fitted edge length
  reaches zero the score has a kink, and R-style finite-difference
  L-BFGS-B stalls before it from every start, as admixtools2's does
  (AADR v66.p1 ancient_5pop admix1: both stop at weight 0.4502-0.4503,
  score 0.67713, against the minimum 0.67706 at 0.45087). Polish the best
  start by default (a bracketing line search on the stalled direction),
  keep the stall under compat.

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
- [x] **H6. qpadm `--compat`** on the CLI, as DEVIATIONS.md step 2 promises.
  AI
  *Done in 6caf8ca:* `popmix qpadm` prints the popdrop table (text,
  JSON) and takes `--compat admixtools2-2.0.10`; `qpgraph --dot` takes it
  for the DOT text. `test_every_compat_switch_is_on_the_cli` holds every
  public function with `compat` to its command.
- [ ] **H7. Calibration by simulation:** a command giving the false-positive
  rate and power of a qpAdm or f4 design on msprime simulations. PR
  Decided: popmix reads simulated genotypes the user makes with msprime
  (GPL-3, never imported); the README's caveats section says so.
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

- [x] **I1. `PopmixError`:** `#[non_exhaustive]`; `PathBuf` paths with
  `#[source]` io errors (`io_with_path`, `error.rs:118`, drops the path except
  for NotFound; bare `?` in `npy.rs`); real variants instead of
  `Validation(String)` (31 uses) and `Numerical(String)` (40); remove the ten
  never-constructed variants; NPY format errors as `InvalidFileFormat`. RS
  [reviewer-verified]
  *Done in 7a6d65e:* `#[non_exhaustive]`; `PathBuf` paths; `Io { path,
  source }` with the io error as `#[source]` and no `From<io::Error>`, so a
  bare `?` cannot drop the path (`IoResultExt::with_path`); `Validation`
  and `Numerical` replaced by `TooFew`, `IndexOutOfRange`,
  `InvalidParameter`, `InvalidData`, `NumericalFailure`, `ThreadPool` and
  `Internal`, population lookups as `PopulationNotFound`; seven
  never-built variants removed; NPY and f2-directory problems are
  `InvalidFileFormat` naming the file. `PopmixError::kind()` (a closed
  `ErrorKind`) is what the bindings map to Python exceptions. The
  seven files that grew moved their tests to `tests.rs` (I7).
- [x] **I2. Typed Python exceptions:** `FileNotFoundError`, `OSError`,
  `KeyError`/`ValueError` for populations and validation, popmix exception
  subclasses; `From<PopmixError> for PyErr` so `.map_err(to_py_err)` becomes
  `?` (everything is `RuntimeError` today, `lib.rs:36`). RS, PD
  *Done in 509a02a:* `popmix.errors`: `PopmixError` and `InputError`,
  `PopulationNotFoundError`, `MissingFileError`, `FileError`,
  `ComputationError`, `InternalError`, each also the built-in it was
  (`ValueError`, ...), so no `except` breaks; not `KeyError`, whose
  `str()` quotes the message. `From<PopmixError> for PyErr` in core behind
  the `python` feature: the 146 `.map_err(to_py_err)` are `?`; the
  bindings' own 33 errors raise the same classes. The CLI keeps the
  traceback of an `InternalError`. README "Errors".
- [x] **I3. Deliberate public API:** private fields with accessors on
  `F2Blocks`/`AlleleFrequencies`; `#[non_exhaustive]` options and results;
  `&Options` not `Option<Options>`; consistent `impl AsRef<Path>`; `impl
  FromStr for GenoFormat`; `rcompat` and the second graph fitter
  `pub(crate)` or documented; distinct names for dates and qpgraph
  `FitOptions`; decide whether ndarray/nalgebra types are public. RS
  *Done in ff23ef1, 5595587 and af93a9e:* `F2Blocks`, `FstBlocks` and
  `AlleleFrequencies` have private fields, accessors and checked
  constructors (`F2Blocks::from_arrays`, `FstBlocks::new`,
  `AlleleFrequencies::new`); the 42 option, compat and result types are
  `#[non_exhaustive]`; the eleven `Option<Options>` parameters are
  `&Options`; popmix-io's path parameters are all `impl AsRef<Path>`
  (prefixes extended with `with_suffix`, never through lossy strings);
  `GenoFormat: FromStr`; the DATES fit options are `DecayFitOptions`.
  The second graph fitter went with C2; `rcompat` stays public, with what
  it promises in its docs. Only ndarray types are public (I5).
- [x] **I4. Move the extraction pipeline out of the bindings:** format
  detection, streaming choice, ploidy detection and relabelling into
  `popmix_io::extract_f2(prefix, &ExtractOptions)` behind a `GenotypeSource`
  trait; `lib.rs` (3,153 lines) becomes thin wrappers, split by domain;
  delegate `best_clade`/`worst_pair` (`lib.rs:2079-2083,2345-2349`) to core.
  RS, AI
  *Done in e389bb0 and e309932:* the bindings' `lib.rs` split by domain
  (3,128 lines to 200; seven modules, each under 800). The pipeline is
  `popmix_io::extract::extract_f2(prefix, &ExtractOptions, warn)`, with
  `popmix_io::dataset` (`open_mapped`, `read_genotypes`, `PloidyRule`,
  `open_vcf`); a function choosing the read path in one place did what a
  `GenotypeSource` trait was meant to, so there is none. Warnings reach
  Python through `popmix_core::Warning` (the start of G13). Found doing it:
  the matrix path ordered populations as the file, the streamed path as
  requested (TGENO and PED gave other `populations` than packed input);
  both follow the request now. `QpWavePairsResult` methods call core's.
- [x] **I5. One linear-algebra library** in public types and hot paths. RS
  *Done in e6d564d:* decided and recorded in popmix-core's crate docs:
  ndarray is the only array type in public signatures (a test keeps it
  so); nalgebra does the small dense solves inside qpAdm, qpGraph and
  ADMIXTURE, faer only PCA's eigendecomposition, measured four times
  nalgebra's speed at 2,000 samples (2.0 s against 8.0 s), the gap
  growing as n^3. No hot loop converts between libraries.
- [x] **I6. Don't promise a stable Rust API at 1.0** until I1-I4 are done
  (Python package 1.0 is separate). RS
  *Done in 2623f61:* README "API stability": semantic versioning
  for the Python package and the CLI from 1.0; the Rust crates stay
  publishable (L1) but promise no stable API, as their crate docs say.
  Open question for the release: crates at 1.0.0 read as stable to Cargo
  users, so publish them as 0.x, or not at all, if they go to crates.io.
- [x] **I7. No source file past about 800 lines,** split by concern, enforced
  by a test. Over it on 2026-10-27: `popmix-python/src/lib.rs` 3,628 (I4),
  `rcompat/lbfgsb.rs` 2,127, `cli/__init__.py` 1,580, `f2_blocks.rs` 1,527,
  `admixture.rs` 1,159, `writers.rs` 1,117, `benchmarks/cases.py` 1,064,
  `dates/ld_decay.rs` 1,003, `streaming.rs` 954, `smartpca.rs` 930,
  `qpgraph_fit.rs` 922, `merge.rs` 848, `qpadm/fit.rs` 845. First move inline
  `mod tests` blocks into sibling `tests.rs` files, then split what is
  still long.
  *Done (2026-10-27):* test modules moved to `tests.rs` (d68dc3d's rule,
  7a6d65e, 97f0764); the bindings' `lib.rs` by domain (e389bb0); L-BFGS-B
  into driver, memory, Cauchy point, subspace and line search (910b747);
  ADMIXTURE into fit, Newton steps and quadratic programs (b9b8d6e); the
  CLI by domain (c0484a7); `benchmarks/cases.py` into a package by area
  (8b58f1a). `test_process_rules.py` keeps every source file under 800.

## J. Python API consistency (P2)

- [x] **J1. One input convention:** `dates()` takes a matrix while `alder`,
  `smartpca`, `king` take a prefix; `qpadm` `sources` vs `qpwave` `left`. PD
  *Done in 2438245:* every option of every public function is
  keyword-only, with its real default in the signature (`fudge=1e-4`, not
  `None`), so 1.x can add options without breaking calls; only the
  populations of a statistic may be given in order, as in admixtools2.
  `dates(prefix, admixed, sources, ...)` reads the populations it needs,
  as every genotype analysis does; `rolloff(prefix, admixed, sources,
  ...)` names its references the same way, and the CLI's `dates` and
  `rolloff` take `--admixed` and `--sources A,B`. The qpAdm family is
  `left, right, target` throughout (admixtools2's `qpadm` order).
  `tests/python/test_api_conventions.py` holds the rules.
- [x] **J2. One result convention:** dicts (`f2`), pyclasses, DataFrames with
  admixtools column names (`est`/`se`); counts as ints, not floats
  (`n_snps: 21077.0`). PD, AI
  Also: `pm.qpwave_pairs` is the Python, admixtools2-shaped table; the
  Rust `_core.qpwave_pairs` and the exported `QpWavePairsResult` are
  returned by no public function. Keep one.
  *Done in dcd788e:* statistics carry admixtools2's names everywhere in
  Python and in the CLI's JSON and TSV: `est, se, z, p, chisq, dof`
  (`weight_se`, `feasible`; qpWave's `left`, `right`; DATES's `admixed`,
  `sources`). Tables are DataFrames; one statistic is a dict with the keys
  of its table row in order (`f2.f3(...)` reads as a row of `qp3pop`),
  now with `p`, counts as ints, z NaN without an SE (null in JSON);
  analyses are result objects. `DatesResult.to_ld_dataframe()` returns
  a DataFrame (it returned a list). The F2Blocks methods `qpdstat`,
  `qp3pop` and `d_stat`, copies of `pm.qpdstat`/`pm.qp3pop` under other
  names, are gone. p-values use musl's erfc, to the last bit of R's.
  Tested in `test_api_conventions.py`.
- [x] **J3. Type hints and NumPy-format docstrings** on the raw PyO3
  functions (CLAUDE.md rule). PD
  *Done in 02c6f2e:* every public function, class and method has a NumPy
  docstring (Parameters, Returns), the Rust bindings' Google-style
  `Args:` included; `popmix/_core.pyi` types the compiled module and
  `py.typed` marks the package typed. `test_docs.py` checks both, the
  stub against the module's real parameter names. Old docs that were
  wrong are corrected (`merge`'s keys, `convertf`'s, `detect_format`'s
  formats, the exception `filter_populations` raises). `extract_f2` and
  `f2_from_frequencies` show `block_size=0.05` rather than None for it.
  `AlleleFrequencies.get` panicked on a SNP index out of range and
  returned None for an unknown population: they raise `IndexError` and
  `PopulationNotFoundError`.
- [x] **J4. pyclass `module="popmix._core"`, names, `frozen` where possible,
  pickling** (reprs show `builtins.PyQpAdmRotateResult`). RS, AI
  *Done in 631ddda:* analysis results are frozen dataclasses in
  `popmix.results` (QpWaveResult, QpAdmResult with PopDropRow,
  F4RatioResult, QpGraphResult with its edges, f3 and starts, ScoredGraph,
  DatesResult with LdBin and DecayFit): they compare by value, pickle, and
  their lists are tuples and arrays read-only, as PcaResult, Sfs,
  RolloffResult and AdmixtureResult already were. The data containers are
  Rust classes named as in Python (`F2Blocks`, `GenotypeMatrix`,
  `AlleleFrequencies`, `Sample`, `Snp`, `AdmixGraph`) in `popmix._core`,
  frozen, and pickle (F2Blocks through its arrays, checked again on
  loading); `Sample`, `Snp` and `GenotypeMatrix` gained constructors.
  `QpAdmResult.is_good_fit` (a p-alone verdict, G1) did not carry over.
  The Rust qpadm_rotate, qpwave_pairs and qpadm_multi bindings and their
  four result classes, used by no public function, are gone (the
  public ones are written in Python).
- [x] **J5. `f4` help text** says "D-statistic / ABBA-BABA" but prints f4. PG
  *Done in 5ef6a52:* the help, the command list and the README say
  f4(pop1, pop2; pop3, pop4) and point to `qpdstat --d-mode` for D
  (`test_cli_help.py`).

## K. Tests, CI and process (P1)

- [ ] **K1. CI test workflow:** `cargo test`, `pytest tests/python`, and the
  golden tests that ship their data (family_chr20_22, rlang, family_4pop) on
  Linux, macOS and Windows, x86_64 and aarch64; a skipped golden fails CI
  (`POPMIX_REQUIRE_GOLDEN=1`). AI, PD
  Decided: Linux for 1.0; macOS and Windows once the repository is public.
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
- [x] **K7. Process rules in CLAUDE.md:** *Done in d68dc3d (`tests/python/test_process_rules.py` checks the mock, `or`-assertion and file-size rules).* a feature is done only when it is
  reachable from the CLI or `popmix.__all__` and an end-to-end test runs it
  on a real result object; no mock result classes when the real type can be
  built; no `or`-chained assertions. AI
- [ ] **K8. Tests that pin the release workflow by running it:** a test tag
  run of `release.yml` on every platform before 1.0.0. AI

## L. Packaging and release (last)

- [x] **L1. popmix-io/core publishable:** `version` on the path dependency,
  `publish = false` on popmix-python, `rust-version` on every crate. RS
  *Done in 5bcaf10:* path dependencies pinned `=1.0.0`, `rust-version` 1.99
  (the pinned toolchain), popmix-python unpublished; `cargo package` of
  both libraries builds. Found doing it: the crate packages lacked LICENSE
  and NOTICE (now symlinked in, packaged as text), and `Cargo.lock` was
  ignored, so wheels were built from freshly resolved versions (now
  tracked, and in the sdist). `tests/python/test_crate_manifests.py`.
- [x] **L2. Unused dependencies:** `ordered-float`; `proptest` and
  `criterion` unless F12/B6 use them. RS *Done in 16dae78:* removed
  `ordered-float` and `proptest` (core), `thiserror` and `approx` (io, found
  by the new test); `criterion` stays for the F12 benches.
  `tests/python/test_dependencies.py` fails on a declared dependency a
  crate does not use.
- [x] **L3. PyO3 upgrade** (0.22 `*_bound` API deprecated from 0.23); abi3
  if rust-numpy allows. RS *Done in f27da1b:* pyo3 0.29.3 and numpy 0.29.0
  (0.30 was under two weeks old; a small step before release); one abi3
  wheel per platform for CPython 3.9 on, checked in release.yml on 3.9 and
  3.14 (`.github/scripts/check-wheel.sh`); `gil_used = true` until popmix is
  tested on free-threaded Python. `convertf`/`merge` moved to `convert.rs`
  and VCF input to `vcf_input.rs` (`lib.rs` 3,628 to 3,152 lines, I7).
- [ ] **L4. Release:** create `trianglegrrl/popmix`, push, run `release.yml`
  by hand, tag `v1.0.0` (PyPI trusted publisher is configured).
- [ ] **L5. Bioconda:** recipe from the PyPI sdist and sha256, local build,
  PR to bioconda-recipes.
- [ ] **L6. nf-core:** pin all 21 modules to the released 1.0.0 container
  (they pin `popmix=0.1.0`, `biocontainers/popmix:0.1.0--py39h2e5f10c_0`);
  add a merge/convertf test from VCF; test data to nf-core/test-datasets;
  modules PR. PG

## M. Minor cleanups (P2)

  Note (J1): the modules call `popmix dates --test-pop/--source1/--source2`
  and `popmix rolloff --refs`; 1.0 takes `--admixed` and `--sources A,B`.
- [x] **M1. Duplicated logic:** three confidence-interval implementations
  (`formatters.py:34`, `rich_output.py:27`, `statistics.py:26`; go with D1).
  AI
  *Done in 53c02fc (D1):* the three modules went with the unused reporting
  layer; the one confidence interval left is READ's (x +- 1.96 SE, as
  READv2), in `kinship/read.rs`. Other copies found since were removed
  where found: the bindings' `best_clade`/`clades` (I4), their population
  reorder for SFS (I3) and their three VCF lookups (16858a9).
- [x] **M2. Idiom:** `f64::midpoint` (9 places), `let … else` instead of
  `unwrap` after `is_none` (`ped.rs:184`, `graph.rs:240`), worthwhile
  `clippy::pedantic` lints (1,782 warnings, mostly noise). RS
  *Done in a44f566:* 21 of clippy's pedantic lints are on for the
  workspace (`[workspace.lints.clippy]`: midpoint, let-else, file
  extensions compared case-blind, `Option<&T>`, needless `Result`,
  unused `self`, hashers, inlined format arguments, `From` for lossless
  casts, ...) and the 288 places they found are fixed; the rest of
  pedantic (casts of counts, doc style, `must_use` on getters) stays off.
  Found doing it: a VCF given with an upper-case extension was not taken
  for one, and a SNP range of a VCF was ignored (fixed before this).
- [x] **M3. Document why graph search is split across Rust and Python** (or
  resolve it with F9). AI
  *Resolved:* it is no longer split. The Python evolutionary search left
  1.0 with C2, so both searches, exhaustive and stepwise, are Rust
  (`qpgraph/graph_search.rs`); `popmix/graph_search.py` only checks the
  method and turns results into `ScoredGraph`s, as its module doc says.
- [x] **M4. CLI consistency** (found writing the skill, N): `sfs` with one
  population prints S/pi/Watterson/Tajima lines to stdout before its TSV
  (even with `-o`), and lacks `--pseudohaploid`; qpadm prints weight SEs
  only in JSON, and qpadm/qpwave have no `--fudge`; `Ignore` is left out by
  `pca`, `kinship`, `admixture` but kept by `extract-f2` without
  `--populations`; `info` orders populations differently for TGENO and
  PACKEDANCESTRYMAP of the same data; `dates` years per generation and
  jackknife, and `rolloff`'s fit range, are fixed in the CLI; JSON from
  `write_table` lacks a trailing newline; README claims `--format json` for
  `kinship` and `sfs` and "6 formats" (TGENO is a seventh), as does the
  convertf docstring.
  *Done in 6b02868:* each as listed, checked by `test_cli_consistency.py`:
  the sfs summary goes to stderr; sfs takes `--pseudohaploid`; sfs and
  kinship take `--format tsv|json`; qpadm prints weight SEs in text and
  tsv; qpadm and qpwave take `--fudge`; dates takes `--generation-time`
  and `--no-jackknife` (its JSON says `admixed`, `sources`); rolloff
  takes `--lval`/`--hval`; extract-f2 leaves out `Ignore`; JSON ends
  with a newline; the docs list seven formats, TGENO with its own README
  section. `info`'s order was already fixed (16858a9). The README's
  qpadm example output turned out invented (M5).
- [x] **M5. README example outputs are real runs.** The qpadm example's
  output was made up, with 3 dof from 3 left and 3 right populations
  (replaced in 6b02868 by a run on AADR v66.p1). Rerun every other
  `**Output:**` block on named data and check them with a test, so they
  cannot drift from what the command prints.
  *Done in 96ee6f3:* every output block is a run on `aadr_v66p1`
  (merge on the merge fixtures, the error example included);
  `tests/golden/test_readme_examples.py` reruns them. It found `info`'s
  population order random (a HashSet in `dataset_info`): first
  appearance now.

## N. Agent skill for popmix (P1)

- [x] **N1. Skill layout with progressive disclosure:** `skills/popmix/SKILL.md`
  (what popmix does, when to use which command, input formats, an index of
  commands) plus one file per CLI command, `skills/popmix/commands/<command>.md`
  (subcommands get their own files), each complete: synopsis, every
  argument and option with type, default and meaning, input and output
  files, worked examples on real data, how to read the result, caveats, and
  the reference tool and function it corresponds to. *Done in d2a252a, at
  `python/popmix/skill/` (package data, the one copy) rather than
  `skills/`.*
- [x] **N2. 1:1 mapping enforced by a test:** walk the Click command tree;
  every command has exactly one doc file and every doc file one command;
  every option and argument of the command appears in its file with the
  right default; fails when the CLI and the skill drift. *Done in d2a252a
  (`tests/python/test_skill.py`).*
- [x] **N3. Examples that run:** the examples in each command file execute in
  the test suite (small bundled data), so they cannot go stale. *Done in
  d2a252a (`# runs` blocks; AADR examples are `# not run: ...`).*
- [x] **N4. Ship it:** include the skill in the package (installed with
  popmix) with a command or documented step to install it for Claude
  Code (`~/.claude/skills/popmix`); version it with popmix. *Done in
  d2a252a (`popmix install-skill`; SKILL.md states the version, tested).*

## O. AADR v66.p1 as the reference dataset (P0)

- [x] **O1. Read TGENO** *Done in 66b3375.* (AdmixTools `transpose_packed`, one row per
  individual) in popmix, streaming; it suits reading a subset of
  individuals. v66.p1 1240K `.geno` begins `TGENO 23089 1233013`.
- [x] **O2. TGENO misdetected as text EIGENSTRAT.** *Done in 66b3375.* `dataset_info` on the v66.p1
  prefix reports `format: eigenstrat`; analyses would read garbage. Detect
  the header and support it (O1) or refuse with a clear error. [confirmed]
- [x] **O3. Build AdmixTools ≥ 8 convertf** *Done in 26e9dbb.* as a reference tool (the system
  EIGENSOFT convertf fails: "bad ID (SNP): TGENO ..."); record its version.
- [x] **O4. Migrate every AADR golden from v62.0 to v66.p1** (`/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB`):
  v66.p1 population labels carry no data-type suffix (`French`, not
  `French.DG`; the suffix is on sample IDs such as `Loschbour.AG`), so
  each reference names its individuals or relabels by suffix. Regenerate
  all references and record the dataset checksums
  (`CHECKSUMS.normalized.md5`).
  *Done (1d771ee..933ebb3 as cherry-picked, then the tolerance commit):
  tests/golden/aadr/build_aadr_v66p1.py builds the labelled dataset from
  checksummed inputs. Four checks the new data broke were resolved with
  measured causes, not looser blanket tolerances: qpgraph_fit stalls at a
  zero-bound kink (G14), near-perfect find_graphs fits get a 1e-12 absolute
  floor, ROLLOFF compares at nls's convergence, and DATES default_result
  skips where DATES-0003 alone exceeds the printed precision.*
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

