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
| **Total** | **257** | **All pass** | **~30 min** | |

Quick tests include:
- DATES validation (7 tests) comparing against DATES v753 reference output
- Merge tests (4 tests) validating SNP intersection and sample combination
- Convertf validation ensuring all 36 format conversions (6x6 matrix) preserve genotypes
- VCF reader tests (11 tests) validating VCF parsing, genotype encoding, and cross-format validation

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
│   ├── test_golden_master.py  # admixtools2 comparison tests
│   └── dates/                 # DATES test data
│       ├── family_chr20_22.*  # Simulated admixture data (54 samples, 21k SNPs)
│       └── chr20_22_reference.json
├── python/                    # Python/CLI tests
│   └── test_cli.py            # CLI command tests
└── rust/                      # Rust integration tests

crates/popmix-io/
└── src/
    ├── merge.rs               # Merge algorithm (12 tests)
    ├── vcf.rs                 # VCF reader (11 tests)
    ├── writers.rs             # Format writers (tested via roundtrip)
    ├── ped.rs                 # PED format reader
    └── ancestrymap.rs         # ANCESTRYMAP format reader
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

## Validation Against Reference Implementations

### DATES vs DATES v753

The `scripts/compare_dates.py` script validates popmix DATES against the original DATES v753:

1. Runs DATES v753 on test data
2. Runs popmix on the same data
3. Compares pair counts (must match **exactly**)
4. Compares correlations (must match within 1e-6)

**Validation results:**
- Total pair count difference: **0** (exact match across all 1000 bins)
- Maximum correlation difference: **4.99e-07** (within 1e-6 tolerance)

### admixtools2 Golden Master Tests

The `tests/golden/test_golden_master.py` tests compare against admixtools2 R package:

| Statistic | Populations | Status |
|-----------|-------------|--------|
| f2 | Han.DG, Yoruba.DG | PASS |
| f2 | Han.DG, Mbuti.DG | PASS |
| FST | Han.DG, Yoruba.DG | PASS |
| FST | Han.DG, Mbuti.DG | PASS |
| FST | Yoruba.DG, Mbuti.DG | PASS |
| f3 | French.DG; Han.DG, Yoruba.DG | PASS |
| f4 | French.DG, Han.DG; Yoruba.DG, Mbuti.DG | PASS |
| qpWave | Multiple populations | PASS |
| f4ratio | Multiple populations | PASS |
| qpAdm | Multiple populations | PASS |

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

**Location**: `tests/golden/*.json`

| File | Reference Implementation | Populations |
|------|-------------------------|-------------|
| `f2_reference.json` | admixtools2 | Han.DG, Yoruba.DG, Mbuti.DG |
| `f3_reference.json` | admixtools2 | French.DG; Han.DG, Yoruba.DG |
| `f4_reference.json` | admixtools2 | French.DG, Han.DG; Yoruba.DG, Mbuti.DG |
| `f4ratio_reference.json` | admixtools2 | Multiple |
| `fst_reference.json` | admixtools2 | Han.DG, Yoruba.DG, Mbuti.DG |
| `qpadm_reference.json` | admixtools2 | Multiple |
| `qpwave_reference.json` | admixtools2 | Multiple |

Generated by `tests/golden/generate_reference.R` using admixtools2 R package.

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
