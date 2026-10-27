# Popmix Development Guide

## MANDATORY: TDD-First Development

**Every feature, bug fix, and refactor starts with a failing test.**

This is not optional. The workflow is:

1. **Write test** that describes expected behavior
2. **Run test** - confirm it fails (RED)
3. **Write code** - minimum to pass
4. **Run test** - confirm it passes (GREEN)
5. **Refactor** - improve while staying green
6. **Commit** - test + implementation together

See `.cursor/rules/008-tdd-first.mdc` for detailed patterns and examples.

**Never commit code without tests. No exceptions.**

## Project Overview

Popmix is a Rust/Python reimplementation of admixtools2 for population genetics
F-statistics and admixture analysis. The project uses:
- **Rust** for core algorithms (popmix-core, popmix-io)
- **PyO3** for Python bindings (popmix-python)
- **Click** for CLI

## Environment Setup (Required)

### Shell Environment

**IMPORTANT**: Conda and the Python venv can conflict. Always unset CONDA_PREFIX before activating the venv:

```bash
# Every session - ALWAYS do this
unset CONDA_PREFIX && source /home/a/popmix/.venv/bin/activate
```

This prevents maturin build errors caused by conflicting Python environments.

### First Time Setup

```bash
cd /home/a/popmix
python -m venv .venv
unset CONDA_PREFIX && source .venv/bin/activate
pip install maturin pytest hypothesis numpy pandas click
```

## Quick Start

```bash
# Build Rust (no venv needed)
cargo build --workspace
cargo test --workspace

# Build + install Python bindings (venv required)
unset CONDA_PREFIX && source .venv/bin/activate
maturin develop -m crates/popmix-python/Cargo.toml

# Run Python tests (venv required)
pytest tests/python -v
```

## Project Structure

```
popmix/
├── crates/
│   ├── popmix-core/     # Pure Rust algorithms
│   ├── popmix-io/       # File format I/O
│   └── popmix-python/   # PyO3 bindings
├── python/popmix/       # Python package
├── tests/
│   ├── rust/            # Rust integration tests
│   ├── python/          # pytest tests
│   └── golden/          # admixtools2 reference data
├── admixtools/          # admixtools2 checkout to run (gitignored); see Provenance
└── .venv/               # Virtual environment (gitignored)
```

## Running R / admixtools2 for Reference

Use the **system R** (not conda R) which has admixtools installed:

```bash
# Start R (system installation, version 4.5)
/usr/bin/R

# Inside R:
library(admixtools)

# Load data
# AADR v66.p1 labelled by data type (Han.DG): build it once with
# python tests/golden/aadr/build_aadr_v66p1.py
prefix <- path.expand("~/.cache/popmix-references/aadr_v66p1/aadr_v66p1")
f2_blocks <- f2_from_geno(prefix, pops = c("Han.DG", "Yoruba.DG", "Mbuti.DG"))

# Compute statistics
fst(f2_blocks)
qpdstat(f2_blocks, pop1 = "Han.DG", pop2 = "Yoruba.DG", pop3 = "Han.DG", pop4 = "Mbuti.DG")
```

### Provenance: never translate GPL source

popmix is MIT-licensed; admixtools2 (GPL-3) and R (GPL-2-or-later) are not.
Implement from papers (Patterson et al. 2012, Harney et al. 2021, Maier et al.
2023, ...), documentation (manuals, `?function` help, vignettes) and the
reference tools' black-box output. Do not translate, or read while writing
code, the source of a GPL tool (admixtools2's R and C++, R's C sources such as
lbfgsb.c and RNG.c, plink-ng): the result would be a derivative work. Running
those tools as oracles is fine, and is how every golden test is made. The
admixtools2 checkout at `/home/a/popmix/admixtools/` (gitignored) is there to
run, and for its documentation (`man/`, vignettes), not to read its code.
Record where a module's method came from in its module doc comment.

## Process Rules

`tests/python/test_process_rules.py` enforces those a machine can check.

- **Done means reachable.** A feature is done only when the CLI or
  `popmix.__all__` reaches it and an end-to-end test runs it there on a real
  result object.
- **No mock result classes** when the real type can be built: tests run the
  real code on small bundled data.
- **One expectation per assertion.** Never `assert a or b`, which passes when
  either half does; assert the exact value, or each part on its own line.
- **Small files, one concern each.** No source file past about 800 lines; put
  new code in a new module rather than growing a long one, and keep long test
  modules in their own files (`#[cfg(test)] mod tests;` with `tests.rs`).
  Files over the limit when the rule came in may only shrink (TASKS.md I7).

## Code Style

### Rust
- No exception swallowing - all errors must propagate with context
- Use `thiserror` for error types
- Prefer `rayon` for parallelization (default: 4 threads)
- Run `cargo fmt` and `cargo clippy` before committing

### Python
- Use type hints everywhere
- Follow PEP 8, enforced by ruff
- Write docstrings in NumPy format

## Testing

### Rust Tests
```bash
cargo test --workspace
cargo test --workspace -- --nocapture  # See output
cargo test -p popmix-core fst          # Run specific tests
```

### Python Tests (always in venv)
```bash
unset CONDA_PREFIX && source .venv/bin/activate
pytest tests/python -v
pytest tests/golden -v      # Compare with admixtools2
```

## Test Data Locations

Real datasets for testing:
- **AADR v66.p1**: `/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB.*`
  (TGENO; labels without data-type suffix, e.g. `French`). The goldens read
  the populations they use, labelled `<label>.<data type>` (`French.DG`), from
  `~/.cache/popmix-references/aadr_v66p1/aadr_v66p1`, which
  `tests/golden/aadr/build_aadr_v66p1.py` builds. Never use AADR v62.0 or
  v66.0: they hold the 161 genomes of Jacobs et al. 2019 that v66.p1 withdrew
  (`tests/python/test_data_ethics.py`).
- **AADR HO**: `/references/aadr_v54.1.p1/HO_public.*` (PLINK format)
- **1000G**: `/8tb/1000genomes/` (PLINK format)

## Numerical Tolerances

When comparing with admixtools2 outputs:
- f-statistics: `rtol=1e-10, atol=1e-15`
- Standard errors: `rtol=1e-6`
- Admixture weights: `rtol=1e-6`
- P-values: `rtol=1e-4`

## Threading

Default: 4 threads. Configurable via:
- Rust: `popmix_core::init_thread_pool(Some(8))`
- Python: `popmix.init_threads(8)`
- CLI: `popmix --threads 8 ...`

## Key Algorithms

### F2-Statistics
```
f2(A,B) = E[(pA - pB)²] - bias_correction
bias_correction = pA(1-pA)/(nA-1) + pB(1-pB)/(nB-1)
```

### F3/F4 (derived from F2)
```
f3(A; B,C) = (f2(A,B) + f2(A,C) - f2(B,C)) / 2
f4(A,B; C,D) = (f2(A,D) + f2(B,C) - f2(A,C) - f2(B,D)) / 2
```

### FST (Fixation Index)
```
FST = mean(f2) / mean(f2_raw + het_a + het_b)
where:
  f2_raw = (pA - pB)²  (no bias correction)
  het = p * (1 - p)    (admixtools2 convention, NOT 2*p*(1-p))
```

Uses ratio-of-means with jackknife standard errors.

## Adding New Features (TDD Workflow)

Follow Test-Driven Development (TDD) with golden master testing:

### 1. RED: Generate Reference Data from admixtools2

First, compute the expected output using admixtools2 in R:

```bash
/usr/bin/R
```

```r
library(admixtools)
# AADR v66.p1 labelled by data type (Han.DG): build it once with
# python tests/golden/aadr/build_aadr_v66p1.py
prefix <- path.expand("~/.cache/popmix-references/aadr_v66p1/aadr_v66p1")
f2 <- f2_from_geno(prefix, pops = c("Han.DG", "Yoruba.DG", "Mbuti.DG"))

# Example: generate FST reference
result <- fst(f2)
print(result)
# Save to JSON for golden master test
```

Save reference data to `tests/golden/<statistic>_reference.json`.

### 2. RED: Write Failing Tests

Write Rust tests that use the reference data:

```rust
#[test]
fn test_fst_matches_admixtools2() {
    let f2_blocks = load_test_f2_blocks();
    let result = f2_blocks.fst("Han.DG", "Yoruba.DG").unwrap();

    // Reference from admixtools2, saved by the generator script
    let expected = load_reference("fst_reference.json")["han_yoruba"].clone();
    // f-statistics: rtol 1e-10, atol 1e-15; standard errors: rtol 1e-6
    assert_relative_eq!(result.estimate, expected.estimate, max_relative = 1e-10, epsilon = 1e-15);
    assert_relative_eq!(result.std_err, expected.std_err, max_relative = 1e-6);
}
```

Run tests to confirm they fail:
```bash
cargo test -p popmix-core fst
```

### 3. GREEN: Implement Minimum Code

Implement the algorithm in Rust from the paper and documentation (see
Provenance above), never from the admixtools2 source.

Run tests until they pass:
```bash
cargo test -p popmix-core fst
```

### 4. Add PyO3 Bindings

Add Python bindings in `crates/popmix-python/src/lib.rs`:

```rust
fn fst(&self, pop_a: &str, pop_b: &str) -> PyResult<Option<HashMap<String, f64>>> {
    // ...
}
```

Rebuild and test:
```bash
unset CONDA_PREFIX && source .venv/bin/activate
maturin develop -m crates/popmix-python/Cargo.toml
pytest tests/python -v
```

### 5. Add Golden Master Test

Add a Python test that compares against saved reference:

```python
def test_fst_golden_master():
    """Compare FST output against admixtools2 reference."""
    f2 = pm.extract_f2(AADR_PREFIX, populations=["Han.DG", "Yoruba.DG"])
    result = f2.fst("Han.DG", "Yoruba.DG")

    with open("tests/golden/fst_reference.json") as f:
        ref = json.load(f)

    est = ref["han_yoruba"]
    np.testing.assert_allclose(result["estimate"], est["estimate"], rtol=1e-10, atol=1e-15)
    np.testing.assert_allclose(result["std_err"], est["std_err"], rtol=1e-6)
```

### 6. REFACTOR

Clean up code while keeping tests green:
- Add error context
- Improve documentation
- Optimize if needed

## F2 Storage Format

F2 blocks are persisted as a directory with:
```
f2_blocks/
├── metadata.json       # Version, populations, block info
├── block_sums.npy      # F2 values (bias-corrected)
├── block_raw_sums.npy  # Raw squared diffs (for FST) [v3+]
├── block_weights.npy   # SNP counts per block
├── block_het_sums.npy  # Heterozygosity sums [v2+]
├── block_het_weights.npy  # Het weights [v2+]
├── block_fst_num.npy   # FST numerator block means, own SNP set/blocks [v4+]
└── block_fst_den.npy   # FST denominator block means [v4+]
```

Current version: 4
