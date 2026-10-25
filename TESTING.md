# Popmix Testing Guide

This document describes how to run tests and validate popmix against reference implementations.

## Quick Start

```bash
# Build and test everything
cargo test --workspace
unset CONDA_PREFIX && source .venv/bin/activate
maturin develop -m crates/popmix-python/Cargo.toml

# Quick tests (includes DATES validation, ~5 min)
pytest tests/ -v

# Full tests (includes AADR integration tests, ~30 min)
pytest tests/ -v --run-slow
```

## Test Summary

| Test Suite | Count | Status | Time | Flag |
|------------|-------|--------|------|------|
| Rust unit tests | 155 | All pass | ~30s | - |
| Rust merge tests | 12 | All pass | ~1s | - |
| Rust convertf tests | 8 | All pass | ~2s | - |
| Rust VCF tests | 11 | All pass | ~1s | - |
| Python quick tests | 26 | All pass | ~3 min | - |
| Python slow tests | 9 | All pass | ~25 min | `--run-slow` |
| Convertf 6x6 matrix | 36 | All pass | ~5s | - |
| Rich output tests | 35 | All pass | ~0.2s | - |
| Formatter tests | 23 | All pass | ~0.5s | - |
| Statistical rigor tests | 24 | All pass | ~0.1s | - |
| HTML report tests | 19 | All pass | ~0.1s | - |
| MultiQC module tests | 16 | All pass | ~0.1s | - |
| Batch mode tests | 17 | All pass | ~0.1s | - |
| **Total** | **391** | **All pass** | **~30 min** | |

Quick tests include:
- DATES validation (7 tests) comparing against DATES v753 reference output
- Merge tests (4 tests) validating SNP intersection and sample combination
- Convertf validation ensuring all 36 format conversions (6x6 matrix) preserve genotypes
- VCF reader tests (11 tests) validating VCF parsing, genotype encoding, and cross-format validation
- Rich output tests (35 tests) for terminal formatting and display
- Formatter tests (23 tests) for JSON, YAML, TSV, Parquet, NDJSON output
- Statistical rigor tests (24 tests) for confidence intervals, significance, effect sizes
- HTML report tests (19 tests) for Plotly-based report generation
- MultiQC module tests (16 tests) for result aggregation
- Batch mode tests (17 tests) for multi-target analysis

## Test Structure

```
tests/
├── conftest.py                # Pytest config (--run-slow option)
├── fixtures/                  # Test fixtures
│   ├── small_eigenstrat.*     # Small dataset (3 samples, 5 SNPs)
│   ├── small_vcf.vcf          # VCF version of small_eigenstrat
│   ├── small_vcf_missing.vcf  # VCF with missing genotypes
│   ├── small_vcf_phased.vcf   # VCF with phased genotypes (|)
│   ├── merge_a.*              # Merge test file A (2 samples, 5 SNPs)
│   └── merge_b.*              # Merge test file B (1 sample, 4 SNPs)
├── golden/                    # Reference data tests
│   ├── test_dates.py          # DATES algorithm tests
│   ├── test_dates_v753.py     # DATES vs DATES v753 runs
│   ├── test_fstats_admixtools2.py  # f-statistics, qpWave, qpAdm vs admixtools2 runs
│   ├── admixtools2/           # admixtools2 references and their generator
│   └── dates/                 # DATES test data
│       ├── family_chr20_22.*  # Simulated admixture data (54 samples, 21k SNPs)
│       └── chr20_22_reference.json
├── python/                    # Python/CLI tests
│   ├── test_cli.py            # CLI command tests
│   ├── test_rich_output.py    # Rich terminal output tests (35 tests)
│   ├── test_formatters.py     # Output format tests (23 tests)
│   ├── test_statistical_rigor.py  # Statistical rigor tests (24 tests)
│   ├── test_html_reports.py   # HTML report generation tests (19 tests)
│   ├── test_multiqc_module.py # MultiQC integration tests (16 tests)
│   └── test_batch_mode.py     # Batch analysis tests (17 tests)
└── rust/                      # Rust integration tests

crates/popmix-io/
└── src/
    ├── merge.rs               # Merge algorithm (12 tests)
    ├── vcf.rs                 # VCF reader (11 tests)
    ├── writers.rs             # Format writers (tested via roundtrip)
    ├── ped.rs                 # PED format reader
    └── ancestrymap.rs         # ANCESTRYMAP format reader

python/popmix/
├── cli/
│   ├── rich_output.py         # Rich terminal formatting
│   ├── formatters.py          # Output format handlers
│   ├── statistics.py          # Statistical rigor utilities
│   └── batch.py               # Batch analysis mode
├── reports/
│   └── html_report.py         # HTML report generation
├── multiqc/
│   ├── module.py              # Main MultiQC module
│   ├── patterns.py            # File pattern detection
│   ├── parsers.py             # JSON file parsing
│   ├── aggregator.py          # Result aggregation
│   └── general_stats.py       # General Statistics columns
└── schemas/
    ├── f2_result.json         # F2 JSON schema
    ├── qpadm_result.json      # qpAdm JSON schema
    └── dates_result.json      # DATES JSON schema
```

## Running Tests

### Rust Tests

```bash
cargo test --workspace
```

Expected output:
```
test result: ok. 111 passed; 0 failed  (popmix-core)
test result: ok. 26 passed; 0 failed   (popmix-io, includes 12 merge tests)
test result: ok. 10 passed; 0 failed   (integration tests)
```

### Merge Tests Specifically

```bash
cargo test -p popmix-io merge
```

Expected output:
```
test merge::tests::test_merge_snp_intersection ... ok
test merge::tests::test_merge_matching_alleles ... ok
test merge::tests::test_merge_flipped_alleles ... ok
test merge::tests::test_strand_ambiguous_detection ... ok
test merge::tests::test_recode_flipped ... ok
test merge::tests::test_merge_combines_samples ... ok
test merge::tests::test_merge_duplicate_sample_error ... ok
test merge::tests::test_merge_preserves_genotypes ... ok
test merge::tests::test_merge_mixed_alleles ... ok
test merge::tests::test_snp_key_hashing ... ok
test merge::tests::test_full_merge_integration ... ok
test merge::tests::test_merge_only_intersection ... ok

test result: ok. 12 passed; 0 failed
```

### VCF Tests Specifically

```bash
cargo test -p popmix-io vcf
```

Expected output:
```
test vcf::tests::test_parse_vcf_gt_hom_ref ... ok
test vcf::tests::test_parse_vcf_gt_het ... ok
test vcf::tests::test_parse_vcf_gt_hom_alt ... ok
test vcf::tests::test_parse_vcf_gt_missing ... ok
test vcf::tests::test_parse_vcf_gt_phased ... ok
test vcf::tests::test_parse_vcf_gt_invalid ... ok
test vcf::tests::test_vcf_reader_basic ... ok
test vcf::tests::test_vcf_reader_missing_data ... ok
test vcf::tests::test_vcf_reader_phased ... ok
test vcf::tests::test_vcf_cross_validation ... ok
test vcf::tests::test_vcf_format_detection ... ok

test result: ok. 11 passed; 0 failed
```

### Python Tests (Quick - includes DATES and Merge)

```bash
unset CONDA_PREFIX && source .venv/bin/activate
pytest tests/ -v
```

Expected output:
```
26 passed in ~3 minutes
```

Includes all DATES validation tests and merge tests.

### Python Tests (Full Suite with AADR)

```bash
pytest tests/ -v --run-slow
```

Expected output:
```
35 passed in ~30 minutes
```

The `--run-slow` flag runs additional integration tests against full AADR reference data (~1.2M SNPs).

### DATES Validation (Scientific-Grade)

```bash
python scripts/compare_dates.py tests/golden/dates/family_chr20_22
```

Expected output:
```
CRITERIA
  Pair counts exact match:     PASS
  Correlations within 1e-6:    PASS

OVERALL: PASS
```

## CLI Testing

### Available Commands

```bash
popmix --help           # Show all commands
popmix info             # Dataset information
popmix populations      # List populations with counts
popmix merge            # Merge two datasets
popmix extract-f2       # Extract F2 blocks
popmix f2               # Compute f2 statistic
popmix f3               # Compute f3 statistic
popmix f4               # Compute f4/D statistic
popmix qpadm            # Admixture modeling
popmix qpwave           # Ancestry wave testing
popmix dates            # Admixture dating
popmix convertf         # Format conversion
```

### Example: Merge Command

```bash
popmix merge \
  --file-a tests/fixtures/merge_a \
  --file-b tests/fixtures/merge_b \
  --output /tmp/merged \
  --output-format eigenstrat
```

Expected output:
```
Merge complete:
  SNPs in file A: 5
  SNPs in file B: 4
  SNPs in common: 3
  Samples in file A: 2
  Samples in file B: 1
  Total samples: 3

Output files written to /tmp/merged
```

### Example: Populations Command

```bash
popmix populations tests/fixtures/small_eigenstrat
```

Expected output:
```
Populations in tests/fixtures/small_eigenstrat (2 total):

Population                               Samples
--------------------------------------------------
PopA                                          2
PopB                                          1

Tip: Use --filter "pattern*" to search
Tip: Specify populations with: --populations "PopA,PopB"
```

### Example: DATES Command

```bash
popmix dates \
  --input tests/golden/dates/family_chr20_22 \
  --test-pop Simulation \
  --source1 French \
  --source2 Yoruba
```

Expected output:
```
DATES Results: Simulation = mix of French and Yoruba
============================================================
Amplitude (fit): 0.036364
Estimated generations: 3.2 ± nan
Estimated years BP: 93 ± nan (29 years/gen)
Fit R²: 0.1634
SNPs used: 20881
LD bins: 1000

LD Decay (first 10 bins):
  Dist(cM)  Pairs      Correlation
     0.1       298384  +0.137219
     0.2       206537  +0.082825
     ...
```

### Example: Population Validation Error

```bash
popmix f2 -i tests/fixtures/small_eigenstrat --pop1 PopA --pop2 WrongPop
```

Expected output:
```
Error: Population 'WrongPop' not found in dataset.

Tip: List all populations with: popmix populations tests/fixtures/small_eigenstrat
Tip: Search populations with: popmix populations tests/fixtures/small_eigenstrat --filter "Wron*"
```

### Example: F2 with AADR Data

```bash
popmix f2 \
  --input /references/aadr/v62.0_1240k_public \
  --pop1 Han.DG \
  --pop2 Yoruba.DG \
  --snp-range 0:10000
```

## Convertf Format Conversion Testing

The `convertf` command converts between 6 genotype file formats. Validation ensures all 36 format pairs (6x6 matrix) work correctly and preserve genotypes.

### Supported Formats

| Format | Extensions | Description |
|--------|------------|-------------|
| `eigenstrat` | `.geno`, `.snp`, `.ind` | Text EIGENSTRAT (0/1/2/9 per SNP) |
| `packedancestrymap` | `.geno`, `.snp`, `.ind` | Binary EIGENSTRAT (2 bits/genotype) |
| `plink` | `.bed`, `.bim`, `.fam` | Binary PLINK |
| `ped` | `.ped`, `.map` | Text PLINK (allele pairs) |
| `ancestrymap` | `.ancestrymapgeno`, `.snp`, `.ind` | Sparse format |
| `vcf` | `.vcf`, `.vcf.gz` | VCF format (plain text or gzip) |

### Quick Validation (6x6 Matrix)

```bash
unset CONDA_PREFIX && source .venv/bin/activate
python3 << 'EOF'
import tempfile
import os
import popmix as pm

prefix = "/home/a/popmix/tests/fixtures/small_eigenstrat"
# Output formats (VCF is input-only for now)
output_formats = ["eigenstrat", "packedancestrymap", "plink", "ped", "ancestrymap"]
# All input formats including VCF
input_formats = output_formats + ["vcf"]

print("6x6 Format Conversion Matrix")
print("=" * 80)

passed = 0
failed = 0
results = []

with tempfile.TemporaryDirectory() as tmpdir:
    # Create source files in each format
    source_files = {}
    for fmt in output_formats:
        out_path = os.path.join(tmpdir, f"source_{fmt}")
        pm.convertf(prefix, out_path, output_format=fmt)
        source_files[fmt] = out_path
    # Add VCF source
    source_files["vcf"] = "/home/a/popmix/tests/fixtures/small_vcf.vcf"

    # Test conversions from all 6 input formats to 5 output formats
    for src_fmt in input_formats:
        for dst_fmt in output_formats:
            src_path = source_files[src_fmt]
            dst_path = os.path.join(tmpdir, f"test_{src_fmt}_to_{dst_fmt}")
            try:
                pm.convertf(src_path, dst_path, input_format=src_fmt, output_format=dst_fmt)
                results.append((src_fmt, dst_fmt, "PASS"))
                passed += 1
            except Exception as e:
                results.append((src_fmt, dst_fmt, f"FAIL: {e}"))
                failed += 1

# Print matrix
header = "FROM / TO"
print(f"\n{header:<20}", end="")
for fmt in output_formats:
    print(f"{fmt[:8]:>10}", end="")
print()
print("-" * 70)

for src_fmt in input_formats:
    print(f"{src_fmt:<20}", end="")
    for dst_fmt in output_formats:
        for r in results:
            if r[0] == src_fmt and r[1] == dst_fmt:
                status = "PASS" if r[2] == "PASS" else "FAIL"
                print(f"{status:>10}", end="")
                break
    print()

print("-" * 70)
print(f"\nResults: {passed}/{passed+failed} conversions passed")
EOF
```

Expected output:
```
6x6 Format Conversion Matrix
================================================================================

FROM / TO             eigenstr  packedan     plink       ped  ancestry
----------------------------------------------------------------------
eigenstrat                PASS      PASS      PASS      PASS      PASS
packedancestrymap         PASS      PASS      PASS      PASS      PASS
plink                     PASS      PASS      PASS      PASS      PASS
ped                       PASS      PASS      PASS      PASS      PASS
ancestrymap               PASS      PASS      PASS      PASS      PASS
vcf                       PASS      PASS      PASS      PASS      PASS
----------------------------------------------------------------------

Results: 30/30 conversions passed
```

Note: VCF is currently supported as an input format only. Output to VCF is not yet implemented.

### CLI Convertf Examples

```bash
# Convert EIGENSTRAT to PLINK
popmix convertf \
  --input tests/fixtures/small_eigenstrat \
  --output /tmp/converted \
  --output-format plink

# Convert with explicit input format
popmix convertf \
  --input /tmp/converted \
  --output /tmp/back_to_eig \
  --input-format plink \
  --output-format eigenstrat

# Verify output files
ls -la /tmp/converted.*
# converted.bed  converted.bim  converted.fam
```

### Format-Specific Notes

**Bit Ordering Differences:**
- PACKEDANCESTRYMAP: Sample 0 in HIGH bits (6-7), encoding: 00=HomRef, 01=Het, 10=HomAlt, 11=Missing
- PLINK: Sample 0 in LOW bits (0-1), encoding: 00=HomAlt, 01=Missing, 10=Het, 11=HomRef

**Genetic Position Units:**
- `.snp` files: Morgans
- `.bim` files: centiMorgans (Morgans * 100)
- `.map` files: centiMorgans

**ANCESTRYMAP Sparse Format:**
- Only stores non-missing genotypes
- One line per genotype: `snp_id sample_id genotype population`

### Rust Unit Tests

```bash
cargo test -p popmix-io writers
cargo test -p popmix-io roundtrip
```

Expected: 8+ writer tests pass

## Merge Algorithm Testing

The merge algorithm is based on PLINK 2.0 by Shaun Purcell & Christopher Chang.

### Test Cases

| Test | Description |
|------|-------------|
| `test_merge_snp_intersection` | Keeps only SNPs present in both files |
| `test_merge_matching_alleles` | Handles identical REF/ALT alleles |
| `test_merge_flipped_alleles` | Recodes genotypes when REF/ALT are swapped |
| `test_strand_ambiguous_detection` | Detects A/T and C/G SNPs |
| `test_merge_combines_samples` | Combines samples from both files |
| `test_merge_duplicate_sample_error` | Errors on duplicate sample IDs |
| `test_merge_preserves_genotypes` | Roundtrip genotype preservation |

### Merge Test Fixtures

**merge_a** (2 samples, 5 SNPs):
```
# merge_a.ind
IND1 M PopA
IND2 F PopA

# merge_a.snp
rs1 1 0.01 100 A G
rs2 1 0.02 200 C T
rs3 1 0.03 300 A C
rs4 2 0.04 400 G T
rs5 2 0.05 500 A G
```

**merge_b** (1 sample, 4 SNPs):
```
# merge_b.ind
IND3 M PopB

# merge_b.snp
rs2 1 0.02 200 C T
rs3 1 0.03 300 A C
rs4 2 0.04 400 G T
rs6 2 0.06 600 A G
```

**Expected merge result**: 3 SNPs (rs2, rs3, rs4), 3 samples (IND1, IND2, IND3)

## Ergonomics and Output Tests

The ergonomics tests validate the user-facing output features:

### Running Ergonomics Tests

```bash
unset CONDA_PREFIX && source .venv/bin/activate
pytest tests/python/test_rich_output.py tests/python/test_formatters.py \
       tests/python/test_statistical_rigor.py tests/python/test_html_reports.py \
       tests/python/test_multiqc_module.py tests/python/test_batch_mode.py -v
```

Expected output:
```
134 passed in ~0.5 seconds
```

### Test Categories

| Test File | Tests | Description |
|-----------|-------|-------------|
| `test_rich_output.py` | 35 | Rich terminal formatting (colors, tables, panels) |
| `test_formatters.py` | 23 | JSON, YAML, TSV, Parquet, NDJSON output |
| `test_statistical_rigor.py` | 24 | CIs, Z-scores, effect sizes, model fit |
| `test_html_reports.py` | 19 | Plotly charts, HTML generation |
| `test_multiqc_module.py` | 16 | File detection, parsing, aggregation |
| `test_batch_mode.py` | 17 | Config parsing, job generation, results |

### Key Test Cases

**Rich Output:**
- Significance coloring (green/yellow/gray based on Z-score)
- Confidence interval display
- qpAdm weight tables with percentages
- DATES timing estimates

**Formatters:**
- JSON schema validation
- YAML human-readability
- Parquet pandas compatibility
- File extension auto-detection

**Statistical Rigor:**
- 95%/99% CI calculation accuracy
- Z-score significance thresholds
- Cohen's d effect size interpretation
- qpAdm weight feasibility assessment

**HTML Reports:**
- Valid HTML5 structure
- Plotly.js integration
- Responsive CSS styling
- qpAdm bar charts, DATES decay curves

**MultiQC:**
- File pattern matching (`*_qpadm.json`, etc.)
- JSON parsing for all result types
- General Statistics column generation
- Sample aggregation

**Batch Mode:**
- YAML config parsing
- Job generation for multiple targets
- Progress tracking
- Summary statistics

## Validation Against Reference Implementations

### DATES vs DATES v753

`tests/golden/test_dates_v753.py` runs popmix on the same data and parameters
as stored DATES v753 runs and checks every stage:

| Check | Compared against | Tolerance |
|-------|------------------|-----------|
| LD table, every bin | `.out` pair counts, covariance, yreg, correlation | pairs exact; others printed precision |
| Full-data fit | `expfit` "mean (generations)" | 5e-4 + DATES simplex tolerance (1e-6 in r = 1e-3 gen) |
| Per-chromosome jackknife | `.jin` | as full-data fit |
| Jackknife estimate and SE | `.jout` | propagated input error, (g + 2) × fit tolerance |
| `default_ld_table` | corrected default vs `.out` | only the bins a correction targets differ; the rest as LD table |
| `default_result` | corrected default vs `.jout` | as jackknife; skipped where a correction changes the LD table |
| `default_vs_exact_mode` | corrected mesh vs DATES exact-mode `.out` | pair counts within 1/qbin per bin |

The first four checks run with `compat="dates-v753"`, which reproduces DATES
exactly. The `default_*` checks cover popmix's corrected default. See
DEVIATIONS.md for what is corrected and why.

`crates/popmix-core/tests/dates_v753_golden.rs` checks the fit and the
weighted jackknife on DATES's own intermediate files, so each stage is
tested on its own.

Configs (`tests/golden/dates/v753/<config>/`):

| Config | Data | qbin | maxdis | DATES v753 | popmix |
|--------|------|------|--------|------------|--------|
| chr20_22_q0 | simulated, chr20-22 | 0 | 1.0 | 1.217 ± 3.648 | matches |
| chr20_22_q10 | simulated, chr20-22 | 10 | 1.0 | 0.751 ± 3.411 | matches |
| simulation_q0 | simulated genome (true: 50) | 0 | 1.0 | 35.680 ± 6.333 | matches |
| simulation_q10 | simulated genome (true: 50) | 10 | 1.0 | 35.738 ± 6.346 | matches |
| aadr_french_q10 | AADR v62: French.DG = Sardinian.DG + Russia_Samara_EBA_Yamnaya.AG | 10 | 0.5 | 140.595 ± 25.037 | matches |

"matches" means `compat="dates-v753"`. The corrected default differs on the
two `_q10` configs with maxdis 1.0 (DATES-0001); DEVIATIONS.md has the numbers.

The simulation and AADR configs are `--run-slow` tests.

DATES fits the correlations it printed to `.out` (6 decimals) and runs the
jackknife on estimates printed to `.jin` (3 decimals). Compat mode does the
same round trip (DATES-0003). On AADR the full-precision default fit differs
by 0.002 generations.

DATES keeps a SNP once each source population has one called sample. It
reads `mincount` from the parameter file but ignores it ("mincount stripped"
in `dates.c`), so popmix defaults `min_source_count` to 1 as well.

**Regenerating references.** Run
`DATES_BIN=<dir> tests/golden/dates/generate_dates_reference.sh [config...]`.
`DATES_BIN` must hold `dates` built with `tests/golden/dates/dates_numbins_fix.patch`,
the DATES helper scripts, and `dowtjack` (upstream builds it only into `src/`).
Unpatched v753 reads an uninitialized `numbins` when sizing buffers, so
`qbin: 0` runs crash or not depending on stack contents. The patch changes
only buffer sizes; mesh-mode output is byte-identical to the unpatched binary.

**Mesh mode (`qbin > 0`) corrects two DATES FFT quirks** (DATES-0001,
DATES-0002): lags past each chromosome's FFT window wrap around, and the
last lag is left unscaled. On chr20-22 with maxdis 1.0 the full-data fit is
3.165 generations in exact mode, 3.160 for popmix's corrected mesh, and
2.890 for DATES's mesh (`compat="dates-v753"`). Exact mode (`qbin = 0`, the
default) has neither quirk.

### Parity and corrections

When a reference tool has a bug that changes results, popmix corrects it by
default and reproduces it under `compat`, so the golden tests can still
require exact agreement. DEVIATIONS.md defines the process (one flag per
deviation; golden, isolation and oracle tests; live numbers) and lists every
deviation.

### f-statistics, qpWave and qpAdm vs admixtools2 2.0.10

`tests/golden/test_fstats_admixtools2.py` (`--run-slow`) runs popmix on the
same AADR v62 populations as stored admixtools2 runs and checks each stage of
the pipeline:

| Check | Compared against | Tolerance |
|-------|------------------|-----------|
| `snp_count` | SNPs left after `maxmiss = 0` and polymorphism filtering | exact |
| `block_lengths` | polymorphic SNPs per 5 cM block | exact, block by block |
| `f2_blocks` | `f2_from_geno` block array, every pair | rtol 1e-10 |
| `f2_estimates` | `f2()` | est rtol 1e-10, se rtol 1e-6 |
| `f3_estimates`, `f4_estimates` | `f3()`, `f4()`, every combination | est rtol 1e-10, se rtol 1e-6 |
| `qpwave_f4`, `qpwave_covariance` | `f2blocks_to_f4stats` | rtol 1e-10, 1e-6 |
| `qpwave_rank_tests` | `qpwave()$rankdrop`, every rank | chisq rtol 1e-6, p rtol 1e-4 |
| `qpadm_weights` | `qpadm()$weights`, unconstrained and constrained | weight and se rtol 1e-6 (atol 1e-10) |
| `qpadm_fit` | `qpadm()$rankdrop` at rank sources - 1 | chisq rtol 1e-6, p rtol 1e-4 |
| `qpadm_popdrop` | `qpadm()$popdrop`, every subset of the sources, popmix run with `compat="admixtools2-2.0.10"` | as `qpadm_weights` and `qpadm_fit`; flags exact |
| `qpadm_popdrop_corrected` | popmix's default nested path (ADMIXTOOLS2-0001): only the nested columns differ from compat, and they follow the nested path walked over admixtools2's own model fits | chisqdiff as chisq, p_nested rtol 1e-4 |
| `f4ratio` | `qpf4ratio()` | alpha and se rtol 1e-6 |
| `qpgraph_f3` | `qpgraph()$f3` est and se, and `ppinv` (inverse f3 covariance) | est rtol 1e-10, se and ppinv rtol 1e-6 |
| `qpgraph_fixed` | `qpgraph` internals at fixed admixture weights (admixtools2's best, and 0.3): edge lengths, fitted f3, score | rtol 1e-10 |
| `qpgraph_fit` | `qpgraph(seed = 1)`, weights optimized | score at most admixtools2's (rtol 1e-10); weights within 1e-6 or admixtools2's spread over its starts |
| `fst_block_lengths`, `fst_blocks`, `fst_estimates` | `f2_from_geno(fst = TRUE)`, `fst()` | exact; rtol 1e-10; est 1e-10, se 1e-6 |

Configs (`tests/golden/admixtools2/v2.0.10/<config>/`):

| Config | Populations | SNPs used | Blocks |
|--------|-------------|-----------|--------|
| dg_6pop | Han.DG, French.DG, Yoruba.DG, Mbuti.DG, Papuan.DG, Chimp.REF | 1,050,751 | 710 |
| ancient_5pop | Russia_Samara_EBA_Yamnaya.AG (pseudo-haploid), Sardinian.DG, French.DG, Mbuti.DG, Chimp.REF | 985,987 | 710 |
| dg_11pop | Druze, Basque, Sardinian, Yoruba, Han, Dinka, Mbuti, Papuan, Karitiana, Ulchi (.DG), Chimp.REF | 1,030,026 | 709 |
| ancient_maxmiss | as ancient_5pop, `maxmiss = 0.25` (one population may be missing) | 1,047,766 | 711 |
| family_4pop | French, Yoruba, Simulation, Simulation.1 (DATES' simulated data, chr 20-22, 54 samples) | 21,080 | 50 |

family_4pop is the only config with fewer than 189 samples, where each
PACKEDANCESTRYMAP record is padded to 48 bytes rather than ceil(n / 4); it
caught popmix's streaming reader assuming the latter.

dg_11pop holds a three-source qpAdm model (Druze = Sardinian + Dinka + Han)
whose constrained solution has two weights inside (0, 1) and one at 0, and a
two-source constrained model inside the bounds, which must equal the
unconstrained one. admixtools2's `qpadm(constrained = TRUE)` fails with three
or more sources (ADMIXTOOLS2-0002), so that model's reference (weights, SEs,
rankdrop, popdrop) comes from the admixtools2 internals `qpadm()` uses (see
DEVIATIONS.md).

qpGraph graphs: a 6-population tree and a 1-admixture graph (dg_6pop), a
1-admixture graph with Yamnaya (ancient_5pop), and a 9-population graph with
two admixtures, one weight on its bound (dg_11pop).

Actual agreement: f2, f3, f4 and FST (1,798 statistics) within 5e-12 relative,
SEs within 1e-12; qpWave chi-squared within 5e-12; qpAdm weights and SEs within
1e-9 (the worst case is an ill-conditioned model with weights 6.06 and -5.06),
model and popdrop chi-squared within 5e-11; f4ratio within 5e-14; qpGraph f3
within 3e-15 (SEs 3e-12), and edge lengths, fitted f3 and score at fixed
weights within 3e-12; with weights optimized from admixtools2's own random
starts, every start's score within 3e-12 and weights within 3e-7, and so the
best weights and their low/high range over the starts.

`crates/popmix-core/tests/rcompat_golden.rs` checks the routines qpGraph
reproduces R with against R 4.5.2 (`tests/golden/rlang/generate_optim.R`):
`runif` after `set.seed` for five seeds (1500 draws each, exact), and
`optim(method = "L-BFGS-B")` on eleven bounded problems, with analytic and
finite-difference gradients, `fnscale`, `lmm`, `factr` and `maxit`: the same
number of evaluations at the same points (1e-8; 1e-4 for a scaled
finite-difference Rosenbrock), the same convergence code and message, and
the result within 1e-10 (1e-8).

`tests/golden/test_find_graphs_admixtools2.py` checks the graph search
against `all_graphs.json` (dg_6pop f2 blocks; 5 populations with 0 and 1
admixtures, 4 with 1 and 2; 2976 graphs):

| Check | Compared against | Tolerance |
|-------|------------------|-----------|
| `test_same_graphs_in_the_same_order` | `generate_all_graphs()`: every graph, in order, by its `graph_hash()` key | exact |
| `test_graph_key_of_admixtools2_graphs` | `graph_hash()` key of each admixtools2 graph | exact |
| `test_qpgraph_at_admixtools2_weights` | `qpgraph()` score of every graph, popmix fitted at admixtools2's best weights | rtol 1e-10 |
| `test_qpgraph_tree_scores` | `qpgraph()` score of every tree (nothing to optimize) | rtol 1e-10 |
| `test_qpgraph_optimized_score_of_every_graph` | `qpgraph()` optimized score of every 1-admixture graph (84 and 1485) | rtol 1e-6 |
| `test_qpgraph_optimized_score_with_two_admixtures` | the same for 2 admixtures: at least 1219 of 1302 agree, at most 33 higher (rounding sends some paths to other minima, DEVIATIONS.md); the best graph's score | rtol 1e-6 |
| `test_exhaustive_search` | `find_graphs(method="exhaustive")` fits the same graphs; its best is at most admixtools2's best | trees rtol 1e-10; best rtol 1e-6 |
| `test_stepwise_search_finds_the_exhaustive_optimum` | stepwise search, 5 populations, 1 admixture | best of all 1485, rtol 1e-6 (an optimized score) |

`tests/golden/test_graph_constraints_admixtools2.py` checks the graph
constraints against `graph_constraints.json` (`generate_graph_constraints.R`):
for every 4th graph of `all_graphs.json` (746 graphs), admixtools2's
`summarize_numadmix`, `summarize_zerof4` and the status (true, false or
undecided) of 6 to 8 event-order constraints of both types, all exactly;
`satisfies_numadmix`, `satisfies_zerof4` and `satisfies_nonzerof4` follow
from them.

`tests/golden/test_graph_search_admixtools2.py` (`--run-slow`) compares the
evolutionary search with admixtools2 `find_graphs()` (`find_graphs.json`,
`generate_find_graphs.R`: 8 seeds each for 5 populations with 1 admixture
and 4 with 2, Chimp.REF as outgroup). Their random numbers differ, so over
30 popmix seeds: the best run reaches the best of admixtools2's runs and the
exhaustive optimum; with one admixture at least 40% of runs reach the
optimum (admixtools2: 5 of 8); and the median number of graphs fitted is
within a factor of 2 of admixtools2's. `tests/python/test_graph_search.py`
checks the graph operations (valid graphs, same leaves, outgroup kept),
reproducibility, constraints, `initgraph` and `max_admix`.

Regenerate with `/usr/bin/Rscript tests/golden/admixtools2/generate_all_graphs.R`
(about 6 minutes).

`tests/golden/test_qpadm_models_admixtools2.py` checks the functions that
evaluate many qpAdm models at once against `qpadm_models.json` (dg_11pop,
`generate_qpadm_models.R`): `qpadm_p` (full-rank models with and without
constraints, a rank-0 model, one without a target), `test_cladality`,
`qpwave_pairs` (5 populations), `qpadm_rotate` and `rotate_models` (25
models), `qpadm_multi` (summary rows and full results) and `qpadm_models`
(15, 16 and 14 models). Models and their order, ranks, degrees of freedom
and feasibility exact; chi-squared rtol 1e-6, p-values 1e-4, weights 1e-6.
Actual agreement: chi-squared within 1.2e-10, p-values 2.9e-8, weights
4.8e-11.

`tests/golden/test_resample_admixtools2.py` checks resampling against
`resample.json` (dg_6pop, `generate_resample.R`): the bootstrap blocks
`boo_list` draws after `set.seed` (exact), `qpgraph_resample_multi` for a
tree and a one-admixture graph and `compare_fits` of their test scores,
`qpgraph_resample_snps` leaving out each of the 710 blocks and with 5
bootstrap resamples, and `qpadm_resample_snps`, `qpwave_resample_snps`,
`qpdstat_resample_snps` and `qp3pop_resample_snps`, replicate by replicate.
Statistics rtol 1e-10, chi-squared and SEs 1e-6, p-values 1e-4, weights 1e-6;
qpGraph scores 1e-8, its weights 1e-6 or within admixtools2's spread over its
starts where the score is flat, out-of-sample scores 1e-6. Actual agreement:
f-statistics within 1.3e-13 (SEs 3e-13), qpAdm weights 1.6e-13, qpGraph
scores 5e-12, out-of-sample scores 4.5e-8, qpGraph weights 1.4e-5 on a graph
whose starts spread over 4e-5. `crates/popmix-core/tests/rcompat_golden.rs`
also checks R's `sample(replace = TRUE)`, which draws the blocks, exactly.

With `--run-slow`, `test_resample_inds` checks the `_resample_inds`
functions against `resample_inds.json` (`generate_resample_inds.R`): three
samples of each dg_6pop population and Chimp.REF; for each of the 15 in a
population of more than one, qpAdm, qpWave, qpdstat, qp3pop and qpGraph on
admixtools2's `f2_from_geno(inds, pops)` without it, at the tolerances
above. `tests/python/test_extract_f2_formats.py` checks `extract_f2(inds,
pops)` itself.

`tests/golden/test_f2_dirs_admixtools2.py` checks admixtools2 f2
directories: popmix reads admixtools2 `extract_f2()`'s directory of dg_6pop
(`f2_dir`, `generate_f2_dir.R`) to admixtools2's f2 blocks exactly, which
match popmix's own (rtol 1e-10, SNP counts exact), and, with `--run-slow`,
admixtools2's `f2_from_precomp()` reads the directory popmix writes to the
same values exactly. `tests/python/test_rds.py` checks the R serialization
underneath.

`tests/golden/test_dot_admixtools2.py` checks `write_dot` against
admixtools2's text (`dot.json`, `generate_dot.R`) for the dg_6pop qpgraph
fits and bare graphs: identical with `compat="admixtools2-2.0.10"`; popmix's
own output, and admixtools2's, parse back to the graph's edges.

`tests/golden/test_classic_parsers.py` checks the parsers of classic
ADMIXTOOLS output against admixtools2's on the output of ADMIXTOOLS 7.0.2
qpDstat, qp3Pop, qpAdm, qpF4ratio, qpGraph and qpfstats
(`tests/golden/admixtools_classic/`: `build_classic.sh`,
`generate_classic.py`; dg_6pop populations of AADR v62): the same tables,
values as read exactly and derived columns within 1e-12.

`tests/golden/test_more_stats_admixtools2.py` checks `lazadm` (constrained
and not), `qpf4diff` (two sets of five populations) and, with `--run-slow`,
`qpfstats` (the fitted f2 blocks of five populations, from AADR genotypes)
against `more_stats.json` (`generate_more_stats.R`): weights and standard
errors rtol 1e-6, f-statistics 1e-10 (qpfstats' blocks 1e-9).

`tests/golden/test_msprime_admixtools2.py` checks that `msprime_sim` and
`msprime_genome` write admixtools2's script exactly for three dg_6pop
graphs with its defaults and other options (`msprime.json`,
`generate_msprime.R`); with `--run-slow` and msprime installed,
`tests/python/test_graph_search.py` runs one.

`tests/golden/test_vcf_ped_admixtools2.py` checks f2 from VCF and PED
input (`vcf_ped/<name>.json`, `generate_vcf_ped.py`: the small family data,
and with `--run-slow` dg_6pop populations of AADR v62, in PLINK format and
exported by PLINK 2 as bgzipped VCF and as PED): popmix's f2 blocks from the
PED (5 cM blocks from its map) and from the VCF (no map: 2 Mb blocks, as
admixtools2 makes without one; populations given by `inds`/`pops`) match
admixtools2 `f2_from_geno()` of the PLINK files at rtol 1e-10, the block
lengths exactly.

### f4-, D- and f3-statistic scans vs admixtools2 2.0.10

`tests/golden/test_dstat_admixtools2.py` checks `qpdstat` and `qp3pop` given
a genotype prefix against admixtools2's `qpdstat()` and `qp3pop()` given the
same prefix (`dstat.json`, from `generate_dstat.R`). These are computed from
the genotypes, not from f2 blocks: every autosomal SNP is used, blocks are
formed over all of them, and each statistic uses the SNPs with data in its
own populations (`allsnps = TRUE`).

| Check | Compared against | Tolerance |
|-------|------------------|-----------|
| `test_ploidy` | ploidy admixtools2 assigns each sample | exact |
| `test_f4_blocks` | `f4blockdat_from_geno(f4mode = FALSE)`: per-block numerator, denominator, SNP count | rtol 1e-10; counts exact |
| `test_f3_blocks` | `f3blockdat_from_geno`: per-block numerator, denominator, estimate, SNP count | rtol 1e-10; counts exact |
| `test_combinations` | the combinations `qpdstat`/`qp3pop` make from a population list, and from vectors | exact, in order |
| `test_qpdstat_f4`, `test_qpdstat_d`, `test_qpdstat_grid_d` | `qpdstat(prefix, ...)`, f4 and D mode | est rtol 1e-10; se and z 1e-6; p 1e-4; n exact |
| `test_qp3pop`, `test_qp3pop_outgroup` | `qp3pop(prefix, ...)`, normalized and `outgroupmode = TRUE` | as qpdstat |
| `test_scans_from_f2_blocks` | `f4()`, `f3()` on f2 blocks (`f3_f4.json`) | est 1e-10, se 1e-6 |

Configs: family_4pop (always run), and with `--run-slow` dg_6pop,
ancient_5pop (a pseudo-haploid population) and mixed_ploidy (AADR samples
relabeled so one population mixes pseudo-haploid and diploid samples, where
the frequencies, which count every sample as diploid, and the f3 correction,
which counts a pseudo-haploid sample once, differ; the test rebuilds the
relabeled .ind beside symlinks to the AADR files). f3 with Chimp.REF (one
pseudo-haploid sample) as target cannot be estimated in either tool.
Actual agreement is within rtol 1e-10 everywhere. Regenerate with
`/usr/bin/Rscript tests/golden/admixtools2/generate_dstat.R [config...]`
(about 95 s per AADR config).

`tests/python/test_dstats.py` checks what is built on the scans, on the small
dataset: combination generation, the outgroup-f3 and D matrices (symmetry,
entries equal to the scans), clustering order, heatmaps, and the `qpdstat`,
`qp3pop`, `outgroup-f3` and `dstat-matrix` commands.

### PCA vs EIGENSOFT 7.2.1 smartpca

`tests/golden/test_pca_eigensoft.py` checks `smartpca` against EIGENSOFT
7.2.1 smartpca (`tests/golden/eigensoft/v7.2.1/<config>/pca.json`). The
stock binary prints eigenvector entries to 4 decimals and eigenvalues to 6,
so the reference comes from smartpca built from the 7.2.1 source with only
its output formats changed to `%.17g` (`smartpca-hiprec.patch`,
`build_smartpca.sh`); `generate_pca.py` runs both and checks that they agree
to the stock binary's printed precision.

| Check | Compared against | Tolerance |
|-------|------------------|-----------|
| `test_snps` | SNPs used in the final pass, SNPs dropped as monomorphic, PCA samples | exact |
| `test_outliers` | `REMOVED outlier` lines: sample, round, eigenvector, score | exact; score rtol 1e-6 |
| `test_eigenvalues` | every eigenvalue (`evaloutname`) | rtol 1e-9 |
| `test_coords` | every sample's coordinates (`evecoutname`), PCA and projected, each eigenvector up to sign | atol 1e-9 |
| `test_tracy_widom` | the Tracy-Widom table, `compat="eigensoft-7.2.1"` | TW and neff rtol 1e-6, p rtol 1e-4 |
| `test_corrected_p_values_differ_only_beyond_the_table` | EIGENSOFT-0001 | exact factor |

Configs: `family_all` (all 54 samples), `family_lsqproject` and
`family_project_no_outliers` (PCA on French and Yoruba, Simulation and
Simulation.1 projected, with and without `lsqproject`), and with
`--run-slow` `aadr_1kg_lsqproject` and `aadr_1kg_project` (PCA on six 1000
Genomes populations, CEU, TSI, FIN, YRI, CHB, GIH; 63 YRI outliers removed
over five rounds; three ancient populations projected). Regenerate with
`tests/golden/eigensoft/build_smartpca.sh`, then
`python tests/golden/eigensoft/generate_pca.py [config...]`.

`tests/python/test_pca.py` checks the API's defaults and errors, the PCA
plot and the `pca` command.

### Kinship vs PLINK 2 (KING-robust) and READv2

`tests/golden/test_king_plink2.py` checks `king` against PLINK v2.0.0-a.7LM
(`tests/golden/plink2/<version>/<config>/king.json`, from
`generate_king.py`): the config's populations are written as PLINK files by
EIGENSOFT convertf, and PLINK 2 run with `--make-king-table counts` (per-pair
counts) and `--make-king bin` (the kinship matrix as doubles).

| Check | Compared against | Tolerance |
|-------|------------------|-----------|
| `test_pairs` | the pairs, in PLINK 2's order | exact |
| `test_counts` | NSNP, HETHET, IBS0 and the het-hom counts (PLINK2-0001: labels swapped) | exact |
| `test_kinship` | `.king.bin`, minus infinity for pseudo-haploid samples | rtol 1e-10 |

Configs: `family` (all 54 samples) and, with `--run-slow`,
`aadr_yri_ceu_lbk` (YRI.DG, CEU.DG and 69 pseudo-haploid Germany_EN_LBK.AG
samples; 36,046 pairs).

`tests/golden/test_read_readv2.py` checks `read_kinship` against READv2 2.01
(`tests/golden/readv2/73739be/<config>/read.json`, from `generate_read.py`,
with READv2 installed as that script describes), run with its defaults on
PLINK files written by EIGENSOFT convertf. With `compat="readv2-2.01"`: the
normalization value (rtol 1e-12), the relationship and first-degree type of
every pair and the overlapping SNPs (exact), and P0, normalized P0, SEs, z
scores, the share of first-degree windows, effective SNPs, the kinship and
its interval (rtol 1e-10). Corrected, P0 and kinship are unchanged and only
the SEs and interval differ. Configs: `family` (49 diploid French and Yoruba
samples) and, with `--run-slow`, `aadr_gurgy` (the 78 pseudo-haploid
Neolithic Gurgy samples, whose pedigree has 59 first-degree pairs: 29 called
parent-offspring, 28 siblings).

`tests/python/test_kinship.py` checks the API's options and the `kinship`
command.

### Site-frequency spectra vs admixtools2 and dadi

`tests/golden/test_sfs.py` checks:

| Check | Compared against | Tolerance |
|-------|------------------|-----------|
| `test_joint_sfs_as_admixtools2` | admixtools2 `joint_sfs()` on `packedancestrymap_to_afs()` (`sfs.json`, from `generate_sfs.R`), `compat="admixtools2-2.0.10"`: every cell of the table, SNPs counted | exact |
| `test_joint_sfs_corrected_keeps_the_inexact_snps` | ADMIXTOOLS2-0003: corrected, every SNP with data in every population is counted | exact |
| `test_projected_sfs_as_dadi` | dadi 2.4.4 `Spectrum.from_data_dict` (`tests/golden/dadi/v2.4.4/<config>/sfs.json`, from `generate_sfs.py`) on allele counts from PLINK 2, projected, unfolded and folded | rtol 1e-10; folded mask exact |
| `test_summary_statistics_as_dadi` | dadi `S`, `pi`, `Watterson_theta`, `Tajima_D` of one-population spectra | rtol 1e-10 |

admixtools2 configs: `family_3pop` (Simulation, Simulation.1, Yoruba) and,
with `--run-slow`, `dg_3pop` (Mbuti.DG, Han.DG, Chimp.REF) and
`ancient_2pop` (pseudo-haploid Yamnaya with missing data, French.DG). dadi
configs: `family` (Yoruba at full size and projected, French, French x
Yoruba) and, with `--run-slow`, `aadr_1kg` (CEU.DG at full size and
projected, YRI.DG, CEU x YRI, CEU x YRI x CHB), autosomes only.

`tests/python/test_sfs.py` checks projection, folding, the summary
statistics' errors, the joint table layout and the `sfs` command.

### ADMIXTURE vs ADMIXTURE 1.3.0

`tests/golden/test_admixture.py` checks `admixture` against ADMIXTURE 1.3.0
(`tests/golden/admixture/v1.3.0/<config>/admixture.json`, from
`generate_admixture.py`), run with `--seed=1 -C=1e-10` on PLINK files
written by EIGENSOFT convertf.

| Check | Compared against | Tolerance |
|-------|------------------|-----------|
| `test_loglik_at_admixture_estimates` | ADMIXTURE's log-likelihood, popmix's evaluated at ADMIXTURE's printed Q (rows renormalized) and P | rtol 1e-9 |
| `test_unsupervised` | final log-likelihood; Q up to cluster order where the maximum is identified | at least ADMIXTURE's; Q atol 1e-4 |
| `test_projection` | `-P` with ADMIXTURE's P fixed (convex) | Q atol 2e-6 |
| `test_supervised` | `--supervised` | at least ADMIXTURE's; Q atol 1e-4 |
| `test_cross_validation` | `--cv=5` (popmix's own folds) | rtol 2% |

Configs: `family_k2`, `family_k3`, `family_projection` (French and Yoruba
fit, Simulation and Simulation.1 projected), `family_supervised`, and with
`--run-slow` `aadr_1kg_k3` (CEU, YRI and CHB, chromosomes 20-22) and
`aadr_1kg_projection` (GIH projected onto it, on the SNPs called in both).

`tests/python/test_admixture.py` checks bounds, the clusters, projection of
fit samples, errors, the bar plot and the `admixture` command.

`tests/python/test_cli.py::TestCommandsMatchAdmixtools2` (`--run-slow`) runs
the CLI itself (`popmix extract-f2`, then `f2`, `f3`, `f4`, `qpwave`, `qpadm`
with `--f2-dir`) on the dg_6pop populations and compares with the same
references at the same tolerances.

`tests/python/test_extract_f2_formats.py` converts the small dataset to
PACKEDANCESTRYMAP, EIGENSTRAT and PLINK with `popmix convertf` and requires
`extract_f2` to give the same f2 blocks from each, and the result dicts and
CLI JSON to have a fixed key order.

`tests/nextflow/test_nfcore_modules.py` (`--run-slow`; needs `nextflow` and
the popmix modules in an nf-core/modules checkout, `NFCORE_MODULES`, default
`~/nf-core/modules`) runs the modules in Nextflow: extract-f2, f2, f3, f4,
qpwave, qpadm and qpgraph on the dg_6pop populations against admixtools2;
findgraphs, dates, convertf and merge on the small dataset against popmix run
directly; and the `geno_fstats_admixture_popmix` subworkflow, whose f3 and f4
are compared with family_4pop. The modules' own nf-tests run with
`nf-test test` in the checkout (24 tests, snapshots included).

FST uses its own SNP set (all SNPs kept by `maxmiss`, monomorphic included)
and blocks, as admixtools2 does, and leaves out Chimp.REF: a single
pseudo-haploid sample has one allele per SNP, so FST is undefined for it
(admixtools2 refuses; popmix returns no FST for its pairs). Regenerate with
`/usr/bin/Rscript tests/golden/admixtools2/generate_reference.R [config...]`.
The SNP set depends on which populations are loaded (`maxmiss = 0` drops a
SNP missing in any of them), so a reference is only valid for its own
population list.

## Test Data

### Merge Test Fixtures

**Location**: `tests/fixtures/merge_a.*` and `tests/fixtures/merge_b.*`

| File | Description |
|------|-------------|
| `merge_a.geno/snp/ind` | 2 samples, 5 SNPs |
| `merge_b.geno/snp/ind` | 1 sample, 4 SNPs (3 overlapping with merge_a) |

**Used by**: Merge CLI tests, merge algorithm validation

### Convertf Test Data (small_eigenstrat)

**Location**: `tests/fixtures/small_eigenstrat.*`

| File | Description |
|------|-------------|
| `small_eigenstrat.geno` | Genotype data (text EIGENSTRAT) |
| `small_eigenstrat.ind` | 3 samples |
| `small_eigenstrat.snp` | 5 SNPs |

**Used by**: Convertf 6x6 matrix validation, roundtrip tests, VCF cross-validation

### VCF Test Data

**Location**: `tests/fixtures/small_vcf*.vcf`

| File | Description |
|------|-------------|
| `small_vcf.vcf` | VCF matching small_eigenstrat (3 samples, 5 SNPs) |
| `small_vcf_missing.vcf` | VCF with missing genotypes (./.) |
| `small_vcf_phased.vcf` | VCF with phased genotypes (0\|1 instead of 0/1) |

**VCF Genotype Encoding**:
- `0/0` → HomRef (0)
- `0/1` or `1/0` → Het (1)
- `1/1` → HomAlt (2)
- `./.` or `.` → Missing (9)
- Phased separator `|` is treated same as unphased `/`

**Used by**: VCF reader tests, convertf matrix validation, cross-format validation

### DATES Test Data (chr20-22)

**Location**: `tests/golden/dates/family_chr20_22.*`

| File | Description |
|------|-------------|
| `family_chr20_22.geno` | Genotype data (PACKEDANCESTRYMAP format) |
| `family_chr20_22.ind` | Sample metadata |
| `family_chr20_22.snp` | SNP positions |

**Contents**:
- 21,162 SNPs across chromosomes 20-22
- 54 samples in 4 populations:
  - `Simulation` (1 sample) - admixed test population
  - `French` - source population 1
  - `Yoruba` - source population 2
  - Reference populations
- Simulated admixture: 20% French, 80% Yoruba, 50 generations ago

**Reference outputs**:
| File | Description |
|------|-------------|
| `chr20_22_reference.json` | DATES v753 output for validation |
| `simulation_reference.json` | Full simulation parameters |
| `Simulation_French_Yoruba/` | Raw DATES v753 output files |

### DATES Reference Output Format

The `chr20_22_reference.json` contains:
```json
{
  "dates_output": {
    "dates_version": 753,
    "coefficient_source1": 0.263
  },
  "ld_decay_bins": [
    {"distance_cm": 0.1, "correlation": 0.137219, "n_pairs": 298384},
    {"distance_cm": 0.2, "correlation": 0.082825, "n_pairs": 206537},
    ...
  ]
}
```

### AADR v62 Reference Data

**Location**: `/references/aadr/v62.0_1240k_public.*`

| File | Size | Description |
|------|------|-------------|
| `.geno` | 5.4 GB | Genotype matrix (PACKEDANCESTRYMAP) |
| `.ind` | 818 KB | 17,629 samples |
| `.snp` | 78 MB | 1,233,013 SNPs |
| `.anno` | 8.2 MB | Sample annotations |

**Format**: PACKEDANCESTRYMAP (EIGENSTRAT binary)

**Used by**: Golden master tests, CLI slow tests

### AADR v54 HO Reference Data

**Location**: `/references/aadr_v54.1.p1/HO_public.*`

| File | Size | Description |
|------|------|-------------|
| `.bed` | 3.0 GB | Genotype matrix (PLINK binary) |
| `.bim` | 28 MB | SNP information |
| `.fam` | 624 KB | Sample information |

**Format**: PLINK binary

**Used by**: CLI PLINK format tests

### Golden Master Reference Files

**Location**: `tests/golden/admixtools2/v2.0.10/<config>/` (generated by
`tests/golden/admixtools2/generate_reference.R`) and
`tests/golden/dates/v753/<config>/` (generated by
`tests/golden/dates/generate_dates_reference.sh`).

| File | Contents |
|------|----------|
| `f2.json` | SNP counts, block lengths, f2 block array, `f2()` |
| `f3_f4.json` | `f3()` and `f4()` for every combination |
| `qpwave.json` | qpWave f4 matrix, covariance, rank tests |
| `qpadm.json` | qpAdm weights, SEs, rankdrop, popdrop (unconstrained and constrained); `qpf4ratio()` |
| `qpgraph.json` | `qpgraph()` fits (edges, score, f3, multistart table, ppinv) and fits at fixed admixture weights |
| `all_graphs.json` (dg_6pop) | every `generate_all_graphs()` graph for four cases: `graph_hash()` key, edges, `qpgraph()` score and best weights |
| `fst.json` | FST block lengths and blocks, `fst()` |
| `qpadm_models.json` (dg_11pop) | `qpadm_p`, `test_cladality`, `qpwave_pairs`, `qpadm_rotate`, `rotate_models`, `qpadm_multi`, `qpadm_models` |
| `resample.json` (dg_6pop) | `boo_list` draws; qpGraph, qpAdm, qpWave, qpdstat and qp3pop `_resample_snps`; `qpgraph_resample_multi`; `compare_fits` |
| `graph_constraints.json` (dg_6pop) | `summarize_numadmix`, `summarize_zerof4` and event-order status for 746 graphs |
| `find_graphs.json` (dg_6pop) | `find_graphs()` best score and graphs fitted, 8 seeds, two cases |
| `more_stats.json` (dg_6pop) | `lazadm`, `qpf4diff`, `qpfstats` |
| `msprime.json` (dg_6pop) | `msprime_sim` and `msprime_genome` scripts for three graphs |
| `resample_inds.json` (dg_6pop) | qpAdm, qpWave, qpdstat, qp3pop and qpGraph with each of 15 samples left out (`f2_from_geno(inds, pops)`) |

### Test Data Requirements by Test Type

| Test Type | Data Required | Location |
|-----------|--------------|----------|
| Merge tests | merge_a, merge_b | `tests/fixtures/` |
| VCF tests | small_vcf*.vcf | `tests/fixtures/` |
| Convertf 6x6 matrix | small_eigenstrat, small_vcf | `tests/fixtures/` |
| DATES validation | chr20-22 data | `tests/golden/dates/` |
| F-statistic golden | AADR v62 | `/references/aadr/` |
| CLI PLINK tests | AADR v54 HO | `/references/aadr_v54.1.p1/` |
| Quick tests | fixtures + dates | `tests/` |
| Slow tests | All above | All locations |

## Numerical Tolerances

| Statistic | Tolerance |
|-----------|-----------|
| Pair counts (DATES) | Exact (0 difference) |
| Correlations (DATES) | 1e-6 |
| F-statistics | rtol=1e-10, atol=1e-15 |
| Standard errors | rtol=1e-6 |
| Admixture weights | rtol=1e-6 |
| P-values | rtol=1e-4 |

## Continuous Integration

```bash
#!/bin/bash
set -e

# Rust tests
cargo test --workspace

# Python setup
python -m venv .venv
source .venv/bin/activate
pip install maturin pytest numpy pandas click

# Build Python bindings
maturin develop -m crates/popmix-python/Cargo.toml

# Python tests (full suite with AADR data)
pytest tests/ -v --run-slow

# DATES validation
python scripts/compare_dates.py tests/golden/dates/family_chr20_22
```

## Debugging Test Failures

### Merge Test Failures

If merge tests fail:
1. Check SNP matching by chromosome + position (not rsID alone)
2. Check allele matching logic (Identical, Flipped, StrandAmbiguous, Mismatch)
3. Check genotype recoding for flipped alleles (0↔2)
4. Verify sample deduplication error handling

### DATES Pair Count Mismatch

If pair counts don't match exactly:
1. Check genetic position storage (should be Morgans, not cM)
2. Check binning formula: `bin = (int)(distance_morgans / 0.001)`
3. Use comparison script to find first diverging pair

### Golden Master Test Failures

1. Check that AADR data exists at `/references/aadr/v62.0_1240k_public`
2. Verify reference JSON files in `tests/golden/`
3. Check population names match exactly (e.g., "Han.DG" not "Han")

### Slow Test Not Running

Make sure to use `--run-slow` flag:
```bash
pytest tests/ -v --run-slow
```

### Convertf Roundtrip Failures

If genotypes differ after roundtrip conversion:

1. **Check bit ordering**: PLINK and PACKEDANCESTRYMAP have opposite bit orders
2. **Check genetic positions**: `.snp` uses Morgans, `.bim`/`.map` use centiMorgans
3. **Check allele encoding**: PED uses ACGT or 1234 encoding
4. **Debug specific conversion**:
   ```python
   import popmix as pm

   # Load original
   orig = pm.read_genotypes("input_prefix")

   # Convert and reload
   pm.convertf("input_prefix", "/tmp/out", output_format="plink")
   conv = pm.read_genotypes("/tmp/out")

   # Compare
   import numpy as np
   orig_arr = orig.to_numpy()
   conv_arr = conv.to_numpy()
   diff_mask = orig_arr != conv_arr
   print(f"Differences: {diff_mask.sum()}")
   print(f"First diff at SNP {np.where(diff_mask.any(axis=1))[0][0]}")
   ```

### Format Detection Issues

If format is not auto-detected correctly:
```python
import popmix as pm

# Check detected format
fmt = pm.detect_format("path/to/prefix")
print(f"Detected: {fmt}")

# Force specific format
pm.convertf("input", "output", input_format="plink", output_format="eigenstrat")
```

### Population Validation Issues

If population validation fails unexpectedly:
```python
import popmix as pm

# List all populations
pop_counts = pm.list_populations("/path/to/data")
for pop, count in pop_counts:
    print(f"  {pop}: {count}")

# Get suggestions for a typo
suggestions = pm.suggest_populations("/path/to/data", "Frnch")
print(suggestions)
```

### VCF Reader Issues

If VCF parsing fails:

1. **Check VCF format version**: Must be VCFv4.x (line starts with `##fileformat=VCFv4`)
2. **Check header line**: Must have `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1...`
3. **Check GT field**: Must be first field in FORMAT column
4. **Check genotype format**: Must be `0/0`, `0/1`, `1/1`, `./.`, or phased equivalents

```python
import popmix as pm

# Test VCF reading
geno = pm.read_vcf("/path/to/file.vcf")
print(f"Samples: {geno.n_samples}")
print(f"SNPs: {geno.n_snps}")

# Check sample names
print(f"Samples: {[s.name for s in geno.samples]}")
```

**Common VCF issues:**
- **Missing FORMAT column**: VCF must have FORMAT field with GT
- **Multi-allelic sites**: Only biallelic SNPs are supported (ALT should be single allele)
- **Gzip detection**: Files ending in `.vcf.gz` are automatically decompressed
- **Sample names**: Taken from header line columns 10+
