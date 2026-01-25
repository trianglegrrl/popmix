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
| Rust unit tests | 135 | All pass | ~30s | - |
| Python quick tests | 20 | All pass | ~5 min | - |
| Python slow tests | 19 | All pass | ~25 min | `--run-slow` |
| **Total Python** | **39** | **All pass** | **~30 min** | |

Quick tests include DATES validation (7 tests) which compares against DATES v753 reference output.

## Test Structure

```
tests/
├── conftest.py                # Pytest config (--run-slow option)
├── golden/                    # Reference data tests
│   ├── test_dates.py         # DATES algorithm tests
│   ├── test_golden_master.py # admixtools2 comparison tests
│   └── dates/                # DATES test data
│       ├── family_chr20_22.* # Simulated admixture data
│       └── chr20_22_reference.json
├── python/                    # Python/CLI tests
│   └── test_cli.py           # CLI command tests
└── rust/                      # Rust integration tests
```

## Running Tests

### Rust Tests

```bash
cargo test --workspace
```

Expected output:
```
test result: ok. 111 passed; 0 failed  (popmix-core)
test result: ok. 14 passed; 0 failed   (popmix-io)
test result: ok. 10 passed; 0 failed   (integration tests)
```

### Python Tests (Quick - includes DATES)

```bash
unset CONDA_PREFIX && source .venv/bin/activate
pytest tests/ -v
```

Expected output:
```
20 passed in ~5 minutes
```

Includes all DATES validation tests using chr20-22 subset data.

### Python Tests (Full Suite with AADR)

```bash
pytest tests/ -v --run-slow
```

Expected output:
```
39 passed in ~30 minutes
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
popmix extract-f2       # Extract F2 blocks
popmix f2               # Compute f2 statistic
popmix f3               # Compute f3 statistic
popmix f4               # Compute f4/D statistic
popmix qpadm            # Admixture modeling
popmix qpwave           # Ancestry wave testing
popmix dates            # Admixture dating
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

### Example: F2 with AADR Data

```bash
popmix f2 \
  --input /references/aadr/v62.0_1240k_public \
  --pop1 Han.DG \
  --pop2 Yoruba.DG \
  --snp-range 0:10000
```

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

### DATES Test Data (chr20-22 subset)

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
| DATES validation | chr20-22 subset | `tests/golden/dates/` |
| F-statistic golden | AADR v62 | `/references/aadr/` |
| CLI PLINK tests | AADR v54 HO | `/references/aadr_v54.1.p1/` |
| Quick tests | chr20-22 only | `tests/golden/dates/` |
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
