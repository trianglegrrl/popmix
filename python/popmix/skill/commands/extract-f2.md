# popmix extract-f2

Compute f2 statistics in jackknife blocks for every pair of the chosen
populations and save them to a directory. Every f-statistic command
(`f2`, `f3`, `f4`, `qpadm`, `qpwave`, `qpgraph`, `find-graphs`, `qpdstat`,
`qp3pop`, `outgroup-f3`, `dstat-matrix`) reads that directory with
`--f2-dir`, so the genotypes are read once per project.

## Synopsis

```text
popmix extract-f2 --input PREFIX --output DIR [--populations A,B,...]
                  [--block-size MORGANS] [--maxmiss FRACTION]
                  [--auto-only | --all-chromosomes] [--pseudohaploid RULE]
                  [--snp-range START:END]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix (any readable format; PACKEDANCESTRYMAP is streamed, others are read into memory first). |
| `--output`, `-o` | text (path) | required | Directory to write; created if needed, files in it replaced. |
| `--populations` | text (comma-separated) | none (all populations in the dataset) | Populations to include. Name them: the blocks hold every pair, so size and time grow with the square of their number, and without the list every population (including `Ignore`) is used. |
| `--block-size` | number | 0.05 | Jackknife block length in Morgans (5 cM), as admixtools2 `blgsize`. Blocks never span chromosomes. |
| `--maxmiss` | number (0 to 1) | 0 | Keep SNPs missing in at most this fraction of the populations (admixtools2 `maxmiss`). 0 keeps only SNPs with data in every population; raise it for ancient populations with low coverage. |
| `--auto-only` / `--all-chromosomes` | flag | `--auto-only` | Use autosomes 1-22 only (default), or every chromosome. |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | Which samples count one allele per SNP. `detect`: samples with no heterozygous call at any SNP; `admixtools2`: none among the first 1,000 SNPs, as admixtools2, which also catches some diploid genomes (ADMIXTOOLS2-0007); `diploid`: none. |
| `--snp-range` | text (`start:end`) | none (all SNPs) | Only SNPs with index start to end-1 (0-based, half-open), read into memory and computed outside the admixtools2 pipeline. Refuses `--maxmiss` and `--pseudohaploid`, and a range holding non-autosomal SNPs unless `--all-chromosomes`. For tests and quick looks, not for analyses. |

## Inputs and outputs

Writes `DIR/` with `metadata.json` (format version 4, populations, number of
blocks, block lengths in SNPs, source prefix) and NumPy arrays:
`block_sums.npy` (bias-corrected f2 per pair and block), `block_raw_sums.npy`,
`block_weights.npy` (SNPs per block), `block_het_sums.npy`,
`block_het_weights.npy`, `block_fst_num.npy`, `block_fst_den.npy` (FST on its
own SNP set). Prints the populations and the number of blocks.

## Examples

```bash
# runs
popmix extract-f2 -i family -o f2_family --populations French,Yoruba,Simulation,Simulation.1
popmix f4 --f2-dir f2_family --pop1 Yoruba --pop2 French --pop3 Simulation --pop4 Simulation.1
popmix qpdstat --f2-dir f2_family --pop1 Yoruba,French,Simulation,Simulation.1
```

```bash
# runs
# every chromosome, a looser missingness filter and 2 cM blocks
popmix extract-f2 -i family6 -o f2_six --all-chromosomes --maxmiss 0.1 --block-size 0.02 \
  --populations French_1,French_2,French_3,Yoruba_1,Yoruba_2,Yoruba_3,Simulation.1
```

```bash
# not run: needs the AADR
popmix --threads 8 extract-f2 -i /references/aadr/v62.0_1240k_public -o f2_project \
  --populations Mbuti.DG,Han.DG,French.DG,Yoruba.DG,Russia_Samara_EBA_Yamnaya,Turkey_N \
  --maxmiss 0.1
```

## Reading the result

`Computed N blocks for K populations`: N is the number of jackknife blocks;
standard errors elsewhere come from leaving out one block at a time. Too few
blocks (under about 20) make standard errors unreliable.

## Caveats

- The SNP set is fixed here: with `--maxmiss 0` it is the SNPs with data in
  all the populations named, and only SNPs polymorphic among them
  (admixtools2 `poly_only = "f2"`). Statistics read from the directory
  later use that set, so results differ slightly from `--input` runs of
  single statistics, which use only the statistic's own populations.
- A block with no data for some pair (possible with `--maxmiss > 0`) is
  dropped for every pair, as admixtools2 `f2_from_geno(remove_na = TRUE)`.
- `--pseudohaploid admixtools2` reproduces admixtools2 exactly
  ([ADMIXTOOLS2-0007](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#admixtools2-0007-diploid-genomes-taken-as-pseudo-haploid));
  the default changes only populations holding a diploid genome whose first
  1,000 calls are homozygous.
- PLINK X, Y, XY and MT chromosomes are read (admixtools2 stops on them,
  [ADMIXTOOLS2-0006](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#admixtools2-0006-plink-chromosome-codes-x-and-y-stop-f2_from_geno)).
- The directory is popmix's format. The Python API reads and writes
  admixtools2 f2 directories (`pm.read_admixtools2_f2`,
  `pm.write_admixtools2_f2`).

## Reference

admixtools2 `extract_f2()` / `f2_from_geno()` (2.0.10); blocks, SNP filters
and values match it to rtol 1e-10.
