# Popmix

A fast Rust/Python toolkit for population genetics analysis. Popmix computes F-statistics, models admixture, and dates population mixing events from ancient and modern DNA data.

## What Does Popmix Do?

Popmix answers questions about human (and other species') population history:

- **How different are two populations?** (f2, FST)
- **Is a population admixed?** (f3, qpAdm, ALDER)
- **Did gene flow occur between populations?** (f4, D-statistics, scans and matrices)
- **What graph of splits and admixtures fits the data?** (qpGraph, graph search)
- **When did populations mix?** (DATES, ALDER, ROLLOFF)
- **What are the ancestry proportions?** (qpAdm, ADMIXTURE)
- **How many ancestral populations contributed?** (qpWave)
- **How are the samples structured and related?** (PCA as smartpca, KING-robust and READ kinship, site-frequency spectra)

Popmix reimplements [admixtools2](https://github.com/uqrmaie1/admixtools), DATES,
smartpca, ADMIXTURE, ALDER and ROLLOFF in Rust with Python bindings and a
command line. Each analysis is checked against the original tool on real data
(TESTING.md). Where an original has a bug, popmix fixes it and offers a
`compat` option that reproduces the original; each one is described in
DEVIATIONS.md.

## Complete Pipeline: VCF to Admixture Analysis

Here's a complete workflow from raw VCF files to population genetics analysis:

```bash
# Step 1: Convert VCF to PLINK format (native VCF support!)
popmix convertf \
  -i my_samples.vcf.gz \
  -o my_samples \
  -of plink

# Step 2: Merge your samples with a reference panel (e.g., AADR)
popmix merge \
  --file-a my_samples \
  --file-b /references/aadr/v62.0_1240k_public \
  --output merged_dataset \
  --output-format plink

# Step 3: See what populations are available
popmix populations merged_dataset
popmix populations merged_dataset --filter "French*"

# Step 4: Extract F2 blocks for the populations you need
popmix extract-f2 \
  -i merged_dataset \
  -o f2_blocks/ \
  --populations "MySample,French.DG,Han.DG,Yoruba.DG,Mbuti.DG"

# Step 5: Run admixture analyses
popmix qpadm \
  -i merged_dataset \
  --target MySample \
  --left "French.DG,Han.DG" \
  --right "Yoruba.DG,Mbuti.DG,Papuan.DG"

# Step 6: Date the admixture event
popmix dates \
  -i merged_dataset \
  --test-pop MySample \
  --source1 French.DG \
  --source2 Han.DG
```

You can also work directly with VCF files for quick analyses:

```bash
# View VCF info directly
popmix info my_samples.vcf.gz

# Merge VCF with reference panel (auto-detects format)
popmix merge \
  --file-a my_samples.vcf.gz \
  --file-b /references/aadr/v62.0_1240k_public \
  --output merged \
  --output-format plink
```

## Quick Start

```bash
# 1. Install
pip install popmix

# 2. Explore your dataset
popmix info /path/to/genotypes
popmix populations /path/to/genotypes

# 3. Compute F2 blocks (the foundation for all statistics)
popmix extract-f2 -i /path/to/genotypes -o f2_blocks/ --populations French,Yoruba,Han

# 4. Run analyses
popmix f2 -i /path/to/genotypes --pop1 French --pop2 Yoruba
popmix qpadm -i /path/to/genotypes --target French --left Yamnaya,Anatolia_N --right Mbuti,Han
```

## Installation

### From PyPI (recommended)

```bash
# Basic installation
pip install popmix

# With CLI enhancements (Rich output, progress bars)
pip install popmix[cli]

# With HTML report generation
pip install popmix[reports]

# With MultiQC integration
pip install popmix[multiqc]

# With plots (matplotlib)
pip install popmix[plot]

```

### From Source

```bash
git clone https://github.com/trianglegrrl/popmix
cd popmix
python -m venv .venv
source .venv/bin/activate
pip install maturin pytest numpy pandas click
maturin develop -m crates/popmix-python/Cargo.toml
```

### From Bioconda

```bash
conda install -c bioconda popmix
```

## Input File Formats

Popmix reads 6 standard genotype file formats. Most formats use a **prefix** convention where related files share a common path prefix. VCF files use a full path (with extension).

### VCF (Variant Call Format)

Standard format from variant callers. Supports both plain text (`.vcf`) and gzip-compressed (`.vcf.gz`) files.

```bash
# VCF uses full path (with extension)
popmix info my_data.vcf.gz
popmix convertf -i my_data.vcf -o my_data_plink -of plink
```

Note: VCF files don't include population information - each sample's population is set to its sample ID by default. You'll typically merge with a reference panel that has population labels.

### EIGENSTRAT (text)

Dense text format with one line per SNP. Files: `.geno`, `.snp`, `.ind`

```
# example.geno - one char per sample (0/1/2/9 for missing)
210012
110022
...

# example.snp - SNP metadata
rs1234  1  0.001  100000  A  G
rs5678  1  0.002  200000  C  T

# example.ind - sample metadata
Sample1  M  French
Sample2  F  Yoruba
```

### PACKEDANCESTRYMAP (binary)

Binary version of EIGENSTRAT. Same `.snp` and `.ind` files, but `.geno` is packed binary with "GENO" header.

### PLINK / PACKEDPED (binary)

Standard PLINK format. Files: `.bed`, `.bim`, `.fam`

```bash
# Widely used, supported by PLINK, ADMIXTOOLS, etc.
data.bed  # Binary genotypes
data.bim  # SNP info (chr, id, cm, bp, a1, a2)
data.fam  # Sample info (fid, iid, pid, mid, sex, pheno)
```

### PED (text)

Text PLINK format. Files: `.ped`, `.map`, `.pedind`

### ANCESTRYMAP (sparse)

Sparse format listing only non-missing genotypes. Files: `.ancestrymapgeno`, `.snp`, `.ind`

```
# example.ancestrymapgeno - one line per genotype
rs1234  Sample1  2
rs1234  Sample2  1
rs5678  Sample1  0
```

### Auto-Detection

Popmix auto-detects the format from file extensions:
- `.vcf` or `.vcf.gz` → VCF (use full path)
- `.bed` present → PLINK
- `.ped` present → PED
- `.ancestrymapgeno` present → ANCESTRYMAP
- `.geno` with "GENO" header → PACKEDANCESTRYMAP
- `.geno` text → EIGENSTRAT

## CLI Commands

### `popmix info`

Display dataset information.

```bash
popmix info /path/to/data
popmix info /path/to/data --populations  # Show all populations with counts
```

**Output:**
```
Dataset: /path/to/data
Format: PACKEDANCESTRYMAP
Samples: 17629
SNPs: 1233013
Populations: 1234
  French.DG: 45 samples
  Han.DG: 32 samples
  Yoruba.DG: 28 samples
  ...
```

### `popmix populations`

List all populations in a dataset with sample counts.

```bash
popmix populations /path/to/data
popmix populations /path/to/data --filter "French*"
popmix populations /path/to/data --format json
```

**Output:**
```
Populations in /path/to/data (847 total):

Population                               Samples
--------------------------------------------------
French.DG                                     45
Han.DG                                        32
Yoruba.DG                                     28
Mbuti.DG                                      12
...

Tip: Use --filter "pattern*" to search
Tip: Specify populations with: --populations "French.DG,Han.DG,Yoruba.DG"
```

### `popmix merge`

Merge two genotype datasets. Keeps only SNPs present in both files.

```bash
popmix merge \
  --file-a my_samples \
  --file-b /references/aadr/v62.0_1240k_public \
  --output merged_data \
  --output-format plink
```

**Output:**
```
Merge complete:
  SNPs in file A: 500000
  SNPs in file B: 1233013
  SNPs in common: 498234
  Samples in file A: 10
  Samples in file B: 17629
  Total samples: 17639
  Flipped SNPs (recoded): 1234

Output files written to merged_data
```

Based on [PLINK 2.0](https://github.com/chrchang/plink-ng) merge algorithm by Shaun Purcell & Christopher Chang.

### `popmix extract-f2`

Extract F2 blocks for downstream analyses. This is the recommended first step for large datasets.

```bash
popmix extract-f2 \
  -i /path/to/genotypes \
  -o f2_blocks/ \
  --populations French,Yoruba,Han,Mbuti \
  --block-size 0.05
```

**Options:**
| Option | Default | Description |
|--------|---------|-------------|
| `--input`, `-i` | required | Input genotype prefix |
| `--output`, `-o` | required | Output directory |
| `--populations` | all | Comma-separated population list |
| `--block-size` | 0.05 | Block size in Morgans |
| `--auto-only` | true | Use only autosomal SNPs |

### `popmix f2`

Compute f2 statistic (genetic drift) between two populations.

```bash
popmix f2 -i /path/to/data --pop1 French --pop2 Yoruba --format json
```

**Output formats:**

Text (default):
```
f2(French, Yoruba):
  Estimate: 0.023456
  Std Error: 0.000123
  Z-score: 190.699
  N SNPs: 1200000
```

JSON (`--format json`):
```json
{
  "estimate": 0.023456,
  "std_err": 0.000123,
  "z_score": 190.699,
  "n_snps": 1200000
}
```

TSV (`--format tsv`):
```
pop1    pop2    estimate    std_err    z_score    n_snps
French  Yoruba  0.023456    0.000123   190.699    1200000
```

### `popmix f3`

Compute f3 statistic (admixture test). A **negative f3** indicates the target is admixed.

```bash
popmix f3 -i /path/to/data --target French --source1 Yamnaya --source2 Anatolia_N
```

**Output:**
```
f3(French; Yamnaya, Anatolia_N):
  Estimate: -0.001234
  Std Error: 0.000089
  Z-score: -13.87
  N blocks: 546
  => Significant admixture signal
```

### `popmix f4`

Compute f4 statistic (D-statistic / ABBA-BABA test). Tests for gene flow.

```bash
popmix f4 -i /path/to/data --pop1 French --pop2 Han --pop3 Yoruba --pop4 Mbuti
```

**Output:**
```
f4(French, Han; Yoruba, Mbuti):
  Estimate: 0.000567
  Std Error: 0.000089
  Z-score: 6.37
  N blocks: 546
  => Significant deviation from tree-like evolution
```

### `popmix qpadm`

Estimate admixture proportions from source populations.

```bash
popmix qpadm \
  -i /path/to/data \
  --target French \
  --left Yamnaya,Anatolia_N,WHG \
  --right Mbuti,Han,Papuan
```

**Output:**
```
qpAdm: French = mix of Yamnaya, Anatolia_N, WHG
Outgroups: Mbuti, Han, Papuan

Admixture weights:
  Yamnaya: 0.4523 (45.2%)
  Anatolia_N: 0.3891 (38.9%)
  WHG: 0.1586 (15.9%)

Model fit: chi-sq=4.2341, df=3, p=0.2371
Feasible: True
=> Good model fit (p > 0.05)
```

### `popmix qpwave`

Test the number of independent ancestry streams.

```bash
popmix qpwave \
  -i /path/to/data \
  --left French,Spanish,Italian \
  --right Mbuti,Han,Papuan
```

**Output:**
```
qpWave: French, Spanish, Italian
Outgroups: Mbuti, Han, Papuan

Rank tests:
  Rank  Chi-sq     df    P-value
     0    45.2341     6  0.000001 *
     1     3.1234     2  0.209876
     2     0.0000     0  1.000000

Estimated rank: 1
=> 2 independent ancestry streams
```

### `popmix dates`

Estimate admixture timing using LD decay.

```bash
popmix dates \
  -i /path/to/data \
  --test-pop French \
  --source1 Yamnaya \
  --source2 Anatolia_N
```

**Output:**
```
DATES Results: French = mix of Yamnaya and Anatolia_N
============================================================
Amplitude (fit): 0.036364
Estimated generations: 150.2 ± 12.3
Estimated years BP: 4356 ± 357 (29 years/gen)
Fit R²: 0.9234
SNPs used: 1200000
LD bins: 1000

LD Decay (first 10 bins):
  Dist(cM)  Pairs      Correlation
     0.1       298384  +0.137219
     0.2       206537  +0.082825
     ...
```

### `popmix convertf`

Convert between genotype file formats.

```bash
# Convert VCF to PLINK (native VCF support!)
popmix convertf -i data.vcf.gz -o data_plink -of plink

# Convert PLINK to EIGENSTRAT
popmix convertf -i data -o data_eigen -of eigenstrat

# Convert with explicit input format
popmix convertf -i data -o data_plink -if eigenstrat -of plink
```

**Input formats:** `vcf`, `ancestrymap`, `eigenstrat`, `packedancestrymap`, `plink`/`packedped`, `ped`

**Output formats:** `ancestrymap`, `eigenstrat`, `packedancestrymap`, `plink`/`packedped`, `ped` (VCF output not supported)

**Format conversion matrix:**

| From \ To | ancestrymap | eigenstrat | packedancestrymap | plink | ped |
|-----------|-------------|------------|-------------------|-------|-----|
| vcf | Yes | Yes | Yes | Yes | Yes |
| eigenstrat | Yes | Yes | Yes | Yes | Yes |
| packedancestrymap | Yes | Yes | Yes | Yes | Yes |
| plink | Yes | Yes | Yes | Yes | Yes |

## Helpful Error Messages

Popmix provides helpful error messages with fuzzy suggestions when population names are wrong:

```
Error: Population 'French' not found in dataset.

Did you mean one of these?
  - French.DG (45 samples)
  - French.SDG (12 samples)
  - French_Alsace.SG (8 samples)

Tip: List all populations with: popmix populations /path/to/data
Tip: Search populations with: popmix populations /path/to/data --filter "French*"
```

## Python API

### Basic Usage

```python
import popmix as pm

# Initialize thread pool (default: 4 threads)
pm.init_threads(8)

# Read genotypes (auto-detects format)
geno = pm.read_genotypes("/path/to/data")
print(f"Loaded {geno.n_snps} SNPs x {geno.n_samples} samples")

# Read VCF directly (use full path with extension)
geno = pm.read_vcf("/path/to/data.vcf.gz")
# Or auto-detect
geno = pm.read_genotypes("/path/to/data.vcf.gz")

# Get population list with counts
pop_counts = pm.list_populations("/path/to/data")
for pop, count in pop_counts[:10]:
    print(f"  {pop}: {count} samples")
```

### Population Discovery

```python
# List all populations with sample counts
pop_counts = pm.list_populations("/path/to/data")

# Get suggestions for a typo
suggestions = pm.suggest_populations("/path/to/data", "Frnch")
# Returns: [("French.DG", 45), ("French.SDG", 12)]
```

### Merging Datasets

```python
# Merge two datasets
result = pm.merge(
    file_a="my_samples",
    file_b="/references/aadr/v62.0_1240k_public",
    output="merged_data",
    output_format="plink"
)
print(f"Merged {result['samples_merged']} samples, {result['snps_merged']} SNPs")
```

### Computing F-statistics

```python
# Method 1: Direct computation (for small datasets)
geno = pm.read_genotypes("/path/to/data")
af = geno.allele_frequencies(populations=["French", "Yoruba", "Han"])
f2 = pm.f2_from_frequencies(af)

# Compute statistics
result = f2.f2("French", "Yoruba")
print(f"f2 = {result['estimate']:.6f} ± {result['std_err']:.6f}")

# Method 2: Streaming extraction (for large datasets, recommended)
f2 = pm.extract_f2("/path/to/data", populations=["French", "Yoruba", "Han"])
f2.save("f2_blocks/")  # Save for later use

# Load precomputed F2 blocks
f2 = pm.read_f2("f2_blocks/")
```

### Admixture Analysis

```python
# qpAdm
result = pm.qpadm(
    f2,
    target="French",
    sources=["Yamnaya", "Anatolia_N"],
    right=["Mbuti", "Han", "Papuan"]
)
print(f"Weights: {dict(zip(result.sources, result.weights))}")
print(f"P-value: {result.p_value}")

# qpWave
result = pm.qpwave(
    f2,
    left=["French", "Spanish", "Italian"],
    right=["Mbuti", "Han", "Papuan"]
)
print(f"Estimated rank: {result.estimated_rank}")
```

### DATES (Admixture Dating)

```python
geno = pm.read_genotypes("/path/to/data")
result = pm.dates(
    geno,
    test_pop="French",
    source1="Yamnaya",
    source2="Anatolia_N"
)
print(f"Admixture: {result.generations:.0f} generations ago")
print(f"Years BP: {result.years_bp:.0f} ± {result.years_bp_se:.0f}")
```

### File Format Conversion

```python
# Convert VCF to PLINK (native VCF support!)
pm.convertf(
    input_prefix="/path/to/data.vcf.gz",  # Full path for VCF
    output_prefix="/path/to/plink_data",
    output_format="plink"
)

# Auto-detect input, specify output
pm.convertf(
    input_prefix="/path/to/plink_data",
    output_prefix="/path/to/eigen_data",
    output_format="eigenstrat"
)

# Explicit formats
pm.convertf(
    input_prefix="/path/to/data",
    output_prefix="/path/to/output",
    input_format="plink",
    output_format="packedancestrymap"
)

# Detect format
fmt = pm.detect_format("/path/to/data.vcf.gz")
print(f"Format: {fmt}")  # "VCF"

fmt = pm.detect_format("/path/to/plink_data")
print(f"Format: {fmt}")  # "PACKEDPED"
```

## Output Formats

All CLI commands support multiple output formats via `--format`:

| Format | Flag | Description |
|--------|------|-------------|
| Text | `--format text` | Human-readable with Rich styling (default) |
| JSON | `--format json` | Machine-readable, schema-validated |
| YAML | `--format yaml` | Human-readable config format |
| TSV | `--format tsv` | Tab-separated, one row per result |
| Parquet | `--format parquet` | High-performance columnar format |
| NDJSON | `--format ndjson` | Newline-delimited JSON for streaming |

### Writing to Files

Use `--output-file` / `-o` to write results directly to files:

```bash
popmix qpadm ... --format json -o results/qpadm_french.json
popmix dates ... --format parquet -o results/dates.parquet
```

The format is auto-detected from the file extension if not specified.

### Rich Terminal Output

All terminal output uses Rich for beautiful, informative displays:
- **Colored significance indicators**: Green (significant), yellow (marginal), gray (not significant)
- **Confidence intervals**: Automatic 95% CI calculation with display
- **Progress bars**: For long-running operations like F2 extraction
- **Formatted tables**: Clean, aligned data presentation

### HTML Reports

Generate standalone, interactive HTML reports with Plotly visualizations:

```bash
popmix qpadm ... --report qpadm_french.html
popmix dates ... --report dates_report.html
```

Reports include:
- Interactive Plotly charts (bar charts, LD decay curves)
- Model fit statistics with visual indicators
- Confidence intervals and uncertainty quantification
- Responsive design for any screen size

## MultiQC Integration

Popmix includes a MultiQC module for aggregating results from multiple analyses:

```bash
# Install with MultiQC support
pip install popmix[multiqc]

# Run MultiQC on popmix output files
multiqc results/ --module popmix
```

The module automatically detects:
- `*_popmix_stats.json` - General statistics
- `*_qpadm.json` - qpAdm results
- `*_dates.json` - DATES results
- `*_f2.json`, `*_f3.json`, `*_f4.json` - F-statistics

MultiQC General Statistics columns include:
- P-value (qpAdm/qpWave model fit)
- Generations and Years BP (DATES)
- R² fit quality
- Number of source populations

## Batch Analysis Mode

Run analyses on multiple targets efficiently:

```bash
# qpAdm batch mode
popmix qpadm-batch \
  -i data \
  --targets targets.txt \
  --left sources.txt \
  --right outgroups.txt \
  --output-dir results/ \
  --report batch_report.html

# DATES batch mode
popmix dates-batch \
  -i data \
  --config dates_config.yaml \
  --output-dir results/
```

Batch features:
- Progress tracking with ETA
- Individual result files per target
- Combined summary JSON/YAML
- Batch summary report with best-fitting models
- Failure tracking and reporting

Example batch config (YAML):
```yaml
targets:
  - French
  - Spanish
  - Italian
sources:
  - Yamnaya
  - Anatolia_N
  - WHG
outgroups:
  - Mbuti
  - Han
  - Papuan
```

## Performance Tips

1. **Use streaming F2 extraction** for large datasets:
   ```bash
   popmix extract-f2 -i data -o f2/ --populations Pop1,Pop2,Pop3
   ```

2. **Specify populations** to avoid computing unnecessary pairs:
   ```python
   f2 = pm.extract_f2(prefix, populations=["A", "B", "C"])  # Fast
   f2 = pm.extract_f2(prefix)  # Slow - computes all pairs
   ```

3. **Use multiple threads**:
   ```bash
   popmix --threads 8 extract-f2 ...
   ```

4. **Use PACKEDANCESTRYMAP format** for fastest I/O

5. **Save F2 blocks** to avoid recomputation:
   ```python
   f2.save("f2_blocks/")
   f2 = pm.read_f2("f2_blocks/")  # Instant reload
   ```

## Validation

Every analysis is checked against its reference tool on real data (AADR v62
and the 1000 Genomes samples in it):

| Analysis | Reference | Agreement |
|----------|-----------|-----------|
| f2, f3, f4, D, FST, scans | admixtools2 2.0.10 | rtol 1e-10 (SEs 1e-6) |
| qpWave, qpAdm, model sets, resampling | admixtools2 2.0.10 | chi-squared and SEs 1e-6, weights 1e-6, p 1e-4 |
| qpGraph | admixtools2 2.0.10 (its random starts and R's L-BFGS-B) | every start's score 3e-12 |
| Graph search | admixtools2 `generate_all_graphs`, `find_graphs` | same graphs and order; optima |
| DATES | DATES v753 | printed precision |
| PCA | EIGENSOFT 7.2.1 smartpca | eigenvalues 1e-9 |
| Kinship | PLINK 2 (KING-robust), READv2 | 1e-10 |
| SFS | admixtools2, dadi 2.4.4 | exact, 1e-10 |
| ADMIXTURE | ADMIXTURE 1.3.0 | log-likelihood at least ADMIXTURE's |
| ALDER | ALDER 1.03 | curves 1e-9, fits 1e-6 |
| ROLLOFF | ADMIXTOOLS 7.0.2 | curves 1e-9; fits at least as good |

See TESTING.md for every check and DEVIATIONS.md for where popmix
corrects the reference tools.

## More analyses

Beyond the commands below, the Python API and CLI cover (each with
`--help` / docstrings):

| Command / function | What |
|--------------------|------|
| `qpdstat`, `qp3pop`, `outgroup-f3`, `dstat-matrix` | f4/D and f3 scans from genotypes or f2 blocks, matrices and heatmaps |
| `qpgraph`, `find-graphs` | admixture graph fits (`--plot`, `--dot`) and searches (stepwise, exhaustive, admixtools2's evolutionary search with its constraints) |
| `pm.qpadm_rotate`, `qpadm_multi`, `qpwave_pairs`, `qpadm_models`, `lazadm` | many qpAdm models at once |
| `pm.*_resample_snps`, `*_resample_inds`, `compare_fits` | jackknife and bootstrap resampling, individuals left out |
| `pca` | smartpca PCA with projection and Tracy-Widom tests |
| `kinship` | KING-robust and READ kinship |
| `sfs` | joint and projected site-frequency spectra |
| `admixture` | ADMIXTURE: unsupervised, supervised, projection, cross-validation |
| `alder`, `rolloff` | admixture dates from weighted LD |
| `pm.read_admixtools2_f2`, `write_admixtools2_f2` | admixtools2 f2 directories |
| `pm.parse_qp*_output` | output of classic ADMIXTOOLS programs |
| `pm.msprime_sim`, `msprime_genome` | msprime simulations of a graph |

## Development

```bash
# Setup
git clone https://github.com/trianglegrrl/popmix
cd popmix
python -m venv .venv
unset CONDA_PREFIX && source .venv/bin/activate
pip install maturin pytest numpy pandas click

# Build
cargo build --workspace
maturin develop -m crates/popmix-python/Cargo.toml

# Test
cargo test --workspace           # Rust tests
pytest tests/ -v                 # Python quick tests (~5 min)
pytest tests/ -v --run-slow      # Full test suite (~30 min)
```

## Citation

If you use Popmix in your research, please cite:

```
Popmix 1.0.0: f-statistics, admixture graphs, qpAdm, PCA, ADMIXTURE and
admixture dating
https://github.com/trianglegrrl/popmix
```

Also cite the original methods:
- F-statistics: Patterson et al. 2012, Genetics
- qpAdm/qpWave: Haak et al. 2015, Nature
- DATES: Narasimhan et al. 2019, Science

## License

MIT License. See LICENSE file for details.
