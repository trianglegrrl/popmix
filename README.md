# Popmix

A fast Rust/Python toolkit for population genetics analysis. Popmix computes F-statistics, models admixture, and dates population mixing events from ancient and modern DNA data.

## What Does Popmix Do?

Popmix answers questions about human (and other species') population history:

- **How different are two populations?** (f2, FST)
- **Is a population admixed?** (f3, qpAdm)
- **Did gene flow occur between populations?** (f4, D-statistics, scans and matrices)
- **What graph of splits and admixtures fits the data?** (qpGraph, graph search)
- **When did populations mix?** (DATES, ROLLOFF)
- **What are the ancestry proportions?** (qpAdm, ADMIXTURE)
- **How many ancestral populations contributed?** (qpWave)
- **How are the samples structured and related?** (PCA as smartpca, KING-robust and READ kinship, site-frequency spectra)

Popmix reimplements [admixtools2](https://github.com/uqrmaie1/admixtools), DATES,
smartpca, ADMIXTURE and ROLLOFF in Rust with Python bindings and a
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
  --file-b /references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB \
  --output merged_dataset \
  --output-format plink

# Step 3: See what populations are available
popmix populations merged_dataset
popmix populations merged_dataset --filter "French*"

# Step 4: Extract F2 blocks for the populations you need
popmix extract-f2 \
  -i merged_dataset \
  -o f2_blocks/ \
  --populations "MySample,French,Han,Yoruba,Mbuti"

# Step 5: Run admixture analyses
popmix qpadm \
  -i merged_dataset \
  --target MySample \
  --left "French,Han" \
  --right "Yoruba,Mbuti,Papuan"

# Step 6: Date the admixture event
popmix dates \
  -i merged_dataset \
  --admixed MySample \
  --sources French,Han
```

You can also work directly with VCF files for quick analyses:

```bash
# View VCF info directly
popmix info my_samples.vcf.gz

# Merge VCF with reference panel (auto-detects format)
popmix merge \
  --file-a my_samples.vcf.gz \
  --file-b /references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB \
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
# Library and command line
pip install popmix

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

### Agent skill for Claude Code

popmix ships a skill for AI agents: `SKILL.md` (which command answers which
question, input formats) and one file per command with every option, its
default, examples and caveats. Install it where Claude Code finds personal
skills, `~/.claude/skills/popmix`:

```bash
popmix install-skill
# after upgrading popmix, replace the old copy
popmix install-skill --force
# or into a project
popmix install-skill --dest .claude/skills/popmix
```

The skill lives in the package (`python/popmix/skill/`); its examples run in
the test suite, so it cannot drift from the CLI.

## Input File Formats

Popmix reads 7 genotype file formats. Most formats use a **prefix** convention where related files share a common path prefix. VCF files use a full path (with extension).

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

### TGENO (binary, by sample)

The packed genotypes of PACKEDANCESTRYMAP stored sample by sample, so that
reading a few samples reads only their records: the format of the AADR from
v66 on (`.geno` with a "TGENO" header, `.snp`, `.ind`).

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

| Command | What |
|---------|------|
| `admixture` | Ancestry proportions with K ancestral populations, as ADMIXTURE. |
| `convertf` | Convert between genotype file formats. |
| `dates` | Estimate admixture date using DATES algorithm. |
| `dstat-matrix` | Matrix of D- (or f4-) statistics for every pair X, Y of populations. |
| `extract-f2` | Extract F2 blocks from genotype data, as admixtools2's extract_f2(). |
| `f2` | Compute f2 statistic between two populations. |
| `f3` | Compute f3 statistic (admixture test). |
| `f4` | Compute the f4 statistic f4(pop1, pop2; pop3, pop4). |
| `find-graphs` | Search for the admixture graphs that best fit the data. |
| `info` | Show dataset information. |
| `install-skill` | Install the popmix agent skill for Claude Code. |
| `kinship` | Kinship for every pair of samples: KING-robust or READ-style. |
| `merge` | Merge two genotype datasets. |
| `outgroup-f3` | Matrix of outgroup f3(outgroup; A, B) for every pair of populations. |
| `pca` | Principal component analysis as EIGENSOFT smartpca, with projection. |
| `populations` | List all populations in a dataset with sample counts. |
| `qp3pop` | f3-statistics f3(pop1; pop2, pop3) for many combinations, as admixtools2 qp3pop. |
| `qpadm` | Estimate admixture weights using qpAdm. |
| `qpdstat` | f4- or D-statistics f(pop1, pop2; pop3, pop4) for many combinations, as admixtools2 qpdstat. |
| `qpgraph` | Fit an admixture graph, as admixtools2 qpgraph does. |
| `qpwave` | Test admixture wave count using qpWave. |
| `rolloff` | ROLLOFF: date admixture from the decay of weighted LD correlation (ADMIXTOOLS). |
| `sfs` | Site-frequency spectrum, as dadi (projected) or admixtools2 joint_sfs(). |

Every command documents its options under `popmix <command> --help`;
the sections below walk through the most used ones. Their outputs are real
runs on `aadr_v66p1`, the 27 populations of AADR v66.p1 the tests use,
labelled by data type (`French.DG`; built by
`tests/golden/aadr/build_aadr_v66p1.py`), and
`tests/golden/test_readme_examples.py` checks that each command still
prints them.

### `popmix info`

Display dataset information.

```bash
popmix info aadr_v66p1
```

**Output:**
```
Dataset: aadr_v66p1
Format: packedancestrymap
Samples: 1201
SNPs: 1233013
Populations: 27
  - Germany_HalberstadtSonntagsfeld_EN_LBK.AG
  - England_MIA.AG
  - Russia_Samara_EBA_Yamnaya.AG
  - Czechia_EBA_CordedWare.AG
  - France_Yonne_N.AG
  ... and 22 more

Tip: Use --populations to show all populations with counts
```

The populations come in the order of the `.ind` file; `popmix info
aadr_v66p1 --populations` lists every one with its sample count.

### `popmix populations`

List all populations in a dataset with sample counts.

```bash
popmix populations aadr_v66p1
```

**Output:**
```
Populations in aadr_v66p1 (27 total):

Population                               Samples
--------------------------------------------------
TSI.DG                                      107
CHB.DG                                      103
YRI.DG                                      102
...
Dinka.DG                                      1

Tip: Use --filter "pattern*" to search
Tip: Specify populations with: --populations "TSI.DG,CHB.DG,YRI.DG"
```

```bash
popmix populations aadr_v66p1 --filter "French*"
```

**Output:**
```
Populations in aadr_v66p1 (1 matched):

Population                               Samples
--------------------------------------------------
French.DG                                    31
```

`--format json` gives the list as JSON.

### `popmix merge`

Merge two genotype datasets. Keeps only SNPs present in both files. Here
on two small EIGENSTRAT files (`tests/fixtures/merge_a` and `merge_b`):

```bash
popmix merge --file-a my_samples --file-b reference --output merged --output-format plink
```

**Output:**
```
File A format (auto-detected): EIGENSTRAT
File B format (auto-detected): EIGENSTRAT
Output format: plink
Merging my_samples + reference -> merged...

Merge complete:
  SNPs in file A: 5
  SNPs in file B: 4
  SNPs in common: 3
  Samples in file A: 2
  Samples in file B: 1
  Total samples: 3

Output files written to merged
```

When alleles had to be swapped, complemented or left out, the summary
counts those SNPs too.

Alleles are reconciled as AdmixTools `mergeit` does: swapped alleles are
recoded, a SNP reported on the other strand is complemented, the unknown
allele of a monomorphic SNP (PLINK `0`, EIGENSTRAT `X`) is filled in, and
A/T, C/G and third-allele SNPs are left out and counted.

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
popmix f2 -i aadr_v66p1 --pop1 French.DG --pop2 Yoruba.DG
```

**Output:**
```
f2(French.DG, Yoruba.DG):
  Estimate: 0.054374
  Std Error: 0.000393
  Z-score: 138.427  (p = 0)
  N SNPs: 1046394
```

```bash
popmix f2 -i aadr_v66p1 --pop1 French.DG --pop2 Yoruba.DG --format json
```

**Output:**
```json
{
  "est": 0.05437365580582565,
  "se": 0.0003927970713996346,
  "z": 138.42683605577471,
  "p": 0.0,
  "n_snps": 1046394
}
```

```bash
popmix f2 -i aadr_v66p1 --pop1 French.DG --pop2 Yoruba.DG --format tsv
```

**Output:**
```
pop1	pop2	est	se	z	p	n_snps
French.DG	Yoruba.DG	0.054374	0.000393	138.427	0	1046394
```

### `popmix f3`

Compute f3 statistic (admixture test). A **negative f3** indicates the target is admixed.

```bash
popmix f3 -i aadr_v66p1 --target French.DG --source1 Russia_Samara_EBA_Yamnaya.AG --source2 Germany_HalberstadtSonntagsfeld_EN_LBK.AG
```

**Output:**
```
f3(French.DG; Russia_Samara_EBA_Yamnaya.AG, Germany_HalberstadtSonntagsfeld_EN_LBK.AG):
  Estimate: -0.001535
  Std Error: 0.000096
  Z-score: -16.064  (p = 4.53e-58)
  N blocks: 713
  => Significant admixture signal
```

### `popmix f4`

Compute the f4 statistic f4(pop1, pop2; pop3, pop4), nonzero when gene flow
makes the two branches share drift. For the D-statistic, use `qpdstat --d-mode`.

```bash
popmix f4 -i aadr_v66p1 --pop1 Mbuti.DG --pop2 French.DG --pop3 Russia_Samara_EBA_Yamnaya.AG --pop4 Germany_HalberstadtSonntagsfeld_EN_LBK.AG
```

**Output:**
```
f4(Mbuti.DG, French.DG; Russia_Samara_EBA_Yamnaya.AG, Germany_HalberstadtSonntagsfeld_EN_LBK.AG):
  Estimate: 0.001193
  Std Error: 0.000188
  Z-score: 6.361  (p = 2e-10)
  N blocks: 713
  => Significant deviation from tree-like evolution
```

### `popmix qpadm`

Estimate admixture proportions from source populations.

```bash
popmix qpadm \
  -i aadr_v66p1 \
  --target French.DG \
  --left Russia_Samara_EBA_Yamnaya.AG,Germany_HalberstadtSonntagsfeld_EN_LBK.AG \
  --right Mbuti.DG,Han.DG,Papuan.DG,Karitiana.DG,Ulchi.DG
```

**Output:**
```
qpAdm: French.DG = mix of Russia_Samara_EBA_Yamnaya.AG, Germany_HalberstadtSonntagsfeld_EN_LBK.AG
Outgroups: Mbuti.DG, Han.DG, Papuan.DG, Karitiana.DG, Ulchi.DG

Admixture weights:
  Russia_Samara_EBA_Yamnaya.AG: 0.4185 (41.8%), SE 0.0275
  Germany_HalberstadtSonntagsfeld_EN_LBK.AG: 0.5815 (58.2%), SE 0.0275

Model fit: chisq=37.7305, dof=3, p=0.0000
Feasible: True
=> Poor model fit (p <= 0.05)

Models (popdrop; pattern 1 = source left out):
  pattern   wt  dof      chisq         p f4rank feasible  best weights
  00         0    3    37.7305    0.0000      1     True  True 0.4185, 0.5815
  01         1    4   434.3187    0.0000      0     True False 1.0000, NA
  10         1    4   196.3585    0.0000      0     True  True NA, 1.0000
```

Two sources are not enough here: present-day French also carry
hunter-gatherer ancestry the model leaves out. The popdrop table fits every
subset of the sources (admixtools2's `popdrop`); `--compat
admixtools2-2.0.10` reproduces admixtools2's nested path (`best`), which
compares siblings rather than nested models (DEVIATIONS.md).

### `popmix qpwave`

Test how many streams of ancestry relate the left populations, as seen from the right ones.

```bash
popmix qpwave \
  -i aadr_v66p1 \
  --left French.DG,Sardinian.DG,Basque.DG \
  --right Mbuti.DG,Han.DG,Papuan.DG,Karitiana.DG,Ulchi.DG
```

**Output:**
```
qpWave: French.DG, Sardinian.DG, Basque.DG
Outgroups: Mbuti.DG, Han.DG, Papuan.DG, Karitiana.DG, Ulchi.DG

Rank tests (* rejected at alpha 0.05):
  Rank   chisq    dof          p
     0  327.8390     8  0.000000 *
     1    5.2071     3  0.157243

Estimated rank: 1 (the lowest rank not rejected)
=> At least 2 streams of ancestry, as seen from these right populations: a
   lower bound, as other right populations or more data can show more streams
   (the rank can be at most 2 here).
```

### `popmix dates`

Estimate admixture timing using LD decay.

```bash
popmix dates \
  -i aadr_v66p1 \
  --admixed French.DG \
  --sources Russia_Samara_EBA_Yamnaya.AG,Germany_HalberstadtSonntagsfeld_EN_LBK.AG \
  --qbin 10
```

**Output:**
```

DATES Results: French.DG = mix of Russia_Samara_EBA_Yamnaya.AG and Germany_HalberstadtSonntagsfeld_EN_LBK.AG
============================================================
Amplitude (fit): 0.003202
Estimated generations: 128.1 ± 21.9
Estimated years BP: 3716 ± 636 (29 years/gen)
Fit R²: 0.3055
SNPs used: 828743
LD bins: 1000

LD Decay (first 10 bins):
  Dist(cM)  Pairs      Correlation
     0.1    718093027  +0.061655
     0.2    612297992  +0.014745
...
```

`--qbin 10` bins SNPs on an FFT mesh, much faster on data this size and
slightly approximate; the default, 0, visits every SNP pair.

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

**Input formats:** `vcf`, `ancestrymap`, `eigenstrat`, `packedancestrymap`, `tgeno`, `plink`/`packedped`, `ped`

**Output formats:** `ancestrymap`, `eigenstrat`, `packedancestrymap`, `tgeno`, `plink`/`packedped`, `ped` (VCF output not supported)

Any input format converts to any output format.

## Helpful Error Messages

Popmix provides helpful error messages with fuzzy suggestions when population names are wrong:

```bash
popmix f2 -i aadr_v66p1 --pop1 French --pop2 Yoruba.DG
```

**Output** (on stderr, exit code 1):
```
Error: Population 'French' not found in dataset.

Did you mean one of these?
  - French.DG (31 samples)

Tip: List all populations with: popmix populations aadr_v66p1
Tip: Search populations with: popmix populations aadr_v66p1 --filter "Fren*"
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
# Returns: [("French", 31), ("FIN", 89), ...]
```

### Merging Datasets

```python
# Merge two datasets
result = pm.merge(
    file_a="my_samples",
    file_b="/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB",
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
print(f"f2 = {result['est']:.6f} ± {result['se']:.6f}")  # also z, p, n_snps

# Method 2: extract_f2 (recommended; PACKEDANCESTRYMAP input is streamed,
# other formats are read into memory first)
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
    left=["Yamnaya", "Anatolia_N"],  # the sources, as admixtools2's left
    right=["Mbuti", "Han", "Papuan"],
    target="French",
)
print(f"Weights: {dict(zip(result.sources, result.weights))}")
print(f"chisq = {result.chisq:.3f}, dof = {result.dof}, p = {result.p:.3g}")

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
result = pm.dates("/path/to/data", "French", ["Yamnaya", "Anatolia_N"])
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

### Errors

Every error popmix raises is a `popmix.PopmixError`, and also the built-in
exception the same failure would be in plain Python, so `except ValueError`
keeps working:

| Exception | Also a | Raised for |
|---|---|---|
| `InputError` | `ValueError` | a bad argument, a malformed file |
| `PopulationNotFoundError` | `ValueError` (an `InputError`) | a population not in the data |
| `MissingFileError` | `FileNotFoundError` | a file that is not there |
| `FileError` | `OSError` | reading or writing a file failed |
| `ComputationError` | `RuntimeError` | valid input, failed computation (a singular system) |
| `InternalError` | `RuntimeError` | a bug in popmix: please report it |

```python
try:
    f2 = pm.extract_f2("data/aadr", ["French", "Frnech"])
except pm.PopulationNotFoundError as e:
    print(e)  # Population not found: 'Frnech'
    print(pm.suggest_populations("data/aadr", "Frnech"))  # [('French', 28), ...]
```

What popmix warns about (VCF records it left out, and later results that
need care) is a `popmix.PopmixWarning`, a `UserWarning`; the command line
prints each as one `Warning: ...` line on stderr. To stop on one instead:

```python
import warnings
warnings.simplefilter("error", pm.PopmixWarning)
```


## Output Formats


Commands that report a few statistics (`f2`, `f3`, `f4`, `qpadm`, `qpwave`,
`dates`, `qpgraph`, ...) print text by default and take `--format json` or
`--format tsv`. The scans (`qpdstat`, `qp3pop`, `outgroup-f3`,
`dstat-matrix`, `kinship`, `sfs`) write a TSV table (or JSON with
`--format json`) to stdout, or to a file with `-o`. `pca`, `admixture`
and `rolloff` write their result files under `--output-prefix`.
Each command's `--help` lists what it takes.

Bad input gives a one-line `Error: ...` and exit code 1; `popmix --debug
...` shows the traceback.

## Performance Tips

1. **Use PACKEDANCESTRYMAP, PLINK .bed or EIGENSTRAT input**: all three are
   memory-mapped and read only for the populations asked for. A VCF is
   decompressed and parsed on every run (packed in memory, 2 bits a
   genotype), so convert a large VCF once:
   ```bash
   popmix convertf -i data -o data_packed -of packedancestrymap
   popmix extract-f2 -i data_packed -o f2/ --populations Pop1,Pop2,Pop3
   ```

2. **Name only the populations you need**: f2 blocks hold every pair of the
   populations given, so their size grows with the square of their number.

3. **Threads** help every analysis, best at about 4 on a laptop: on a
   machine with few performance cores, 16 threads can be slower than 4:
   ```bash
   popmix --threads 4 pca ...
   ```

4. **Save F2 blocks** to avoid recomputation:
   ```python
   f2.save("f2_blocks/")
   f2 = pm.read_f2("f2_blocks/")  # Instant reload
   ```

## Benchmarks

popmix against each reference tool on the same 473 AADR v66.p1 individuals
(1.23M SNPs; PCA and kinship on their 220 present-day samples), median of 5
runs after a warm-up, on a 16-thread laptop (i7-1360P, 62 GB):

| Analysis | Reference | Reference | popmix | Agreement |
|----------|-----------|-----------|--------|-----------|
| convertf (packed to text) | AdmixTools 8 convertf | 78.0 s | 1.6 s | identical files |
| f2 blocks, 16 populations | admixtools2 `f2_from_geno` | 36.9 s, 7.9 GB | 0.9 s, 0.7 GB | every entry 1e-10 |
| qpAdm from genotypes | admixtools2 `qpadm` | 21.5 s | 1.0 s | 2.0e-11 |
| qpAdm rotation, 30 models | admixtools2 `qpadm_rotate` | 21.4 s | 1.3 s | 6.9e-10 |
| qpDstat, 20 f4 | admixtools2 `qpdstat` | 16.8 s | 1.6 s | 3.4e-11 |
| qpGraph, 10 starts | admixtools2 `qpgraph` | 16.8 s | 0.9 s | 2.8e-12 |
| Graph search | admixtools2 `find_graphs` | 45.2 s | 2.2 s | same optimum |
| qpGraph bootstrap | admixtools2 `qpgraph_resample_snps` | 18.9 s | 1.3 s | 4.1e-12 |
| qpDstat, qp3Pop, qpGraph | ADMIXTOOLS 7.0.2 | 42.9 s, 2.4 min, 91.0 s | 1.5 s, 1.5 s, 0.8 s | printed precision |
| PCA, 10 PCs | EIGENSOFT smartpca | 83.2 s | 4.3 s | 4.3e-12 |
| ROLLOFF | ADMIXTOOLS 7.0.2 rolloff | 4.3 min | 15.7 s | 4.1e-12 |
| DATES | DATES v753 | 30.7 s | 2.3 s | printed precision |
| READ kinship | READv2 | 28.3 s | 0.8 s | 4.8e-12 |
| KING-robust kinship | PLINK 2 (16 threads) | 1.2 s | 1.1 s | 5.0e-13 |
| ADMIXTURE K=3 | ADMIXTURE 1.3.0 | 13.9 s | 3.6 s | log-likelihood 0.004 higher |

popmix is faster than the reference in every case (KING only just: at
3,970 samples it is 7% behind PLINK 2, TASKS.md F14). On one thread
popmix's ADMIXTURE K=3 fit takes 10.0 s, ADMIXTURE 14.2 s. Its output is byte-identical at 1, 4 and 16 threads in every case. Classic qpAdm
differs from ADMIXTOOLS 7.0.2 by design (DEVIATIONS.md).

Scaling (`benchmarks/scaling.py`): f2 blocks of 497 populations take 5.8 s
and 8.2 GB, and of the 16 populations from all 23,089 AADR samples 1.0 s
(admixtools2: 67 s). At 3,970 samples PCA takes 8.6 min (smartpca 23 min)
and KING 3.6 min (PLINK 2 3.4 min). Full results, with CPU time, memory and
run-to-run variation, are in `benchmarks/results/`.

## Validation

Every analysis is checked against its reference tool on real data (AADR
v66.p1 and the 1000 Genomes samples in it):

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
| ROLLOFF | ADMIXTOOLS 7.0.2 | curves 1e-9; fits at least as good |

See TESTING.md for every check and DEVIATIONS.md for where popmix
corrects the reference tools.

## More analyses

Beyond the walkthroughs above, the Python API and CLI cover (each with
`--help` / docstrings):

| Command / function | What |
|--------------------|------|
| `qpdstat`, `qp3pop`, `outgroup-f3`, `dstat-matrix` | f4/D and f3 scans from genotypes or f2 blocks, matrices and heatmaps |
| `qpgraph`, `find-graphs` | admixture graph fits (`--plot`, `--dot`) and searches (stepwise, exhaustive) |
| `pm.qpadm_rotate`, `qpadm_multi`, `qpwave_pairs`, `qpadm_models` | many qpAdm models at once |
| `pm.*_resample_snps`, `*_resample_inds`, `compare_fits` | jackknife and bootstrap resampling, individuals left out |
| `pca` | smartpca PCA with projection and Tracy-Widom tests |
| `kinship` | KING-robust and READ kinship |
| `sfs` | joint and projected site-frequency spectra |
| `admixture` | ADMIXTURE: unsupervised, supervised, projection, cross-validation |
| `rolloff` | admixture dates from weighted LD |
| `pm.read_admixtools2_f2`, `write_admixtools2_f2` | admixtools2 f2 directories |

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

## API stability

From 1.0, popmix follows semantic versioning for its Python package (the
names in `popmix.__all__`, their parameters and their results) and its
CLI (commands, options and output formats): a 1.x release adds to them
but does not break them. Analyses may still correct a result against its
reference tool; each such change is listed in DEVIATIONS.md.

The Rust crates (`popmix-core`, `popmix-io`, `popmix-python`) are
popmix's internals. The two libraries can be built on, but their Rust API
has no stability promise: it may change in any release, 1.x included.

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
