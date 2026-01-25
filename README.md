# Popmix

A Rust/Python reimplementation of admixtools2 for population genetics F-statistics and admixture analysis.

## Features

- **F-statistics**: f2, f3, f4, Fst computation
- **Admixture analysis**: qpAdm, qpWave, qpGraph
- **File formats**: PLINK, EIGENSTRAT, PACKEDANCESTRYMAP
- **Performance**: Rust backend with parallel computation (default 4 threads)
- **Python API**: Native Python bindings via PyO3

## Installation

```bash
# From source
pip install maturin
maturin develop

# Or from bioconda (when available)
conda install -c bioconda popmix
```

## Quick Start

```python
import popmix as pm

# Initialize with 4 threads (default)
pm.init_threads(4)

# Load data and compute f-statistics
# (Full API coming soon)
```

## CLI Usage

```bash
popmix extract-f2 --input data/genotypes --output f2_dir/
popmix f4 --data f2_dir/ --pop1 A --pop2 B --pop3 C --pop4 D
popmix qpadm --data f2_dir/ --target French --left "WHG,EEF" --right "Mbuti,Chimp"
```

## Development

```bash
# Set up venv (required)
python -m venv .venv
source .venv/bin/activate
pip install maturin pytest numpy pandas

# Build
cargo build --workspace
maturin develop -m crates/popmix-python/Cargo.toml

# Test
cargo test --workspace
pytest tests/
```

## License

MIT
