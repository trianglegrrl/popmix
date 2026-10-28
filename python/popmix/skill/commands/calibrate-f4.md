# popmix calibrate-f4

Compute one f4 statistic on each of many simulated datasets and count how often it is significantly non-zero: on a simulated tree, its false-positive rate; with simulated gene flow, its power.

## Synopsis

```text
popmix calibrate-f4 REPLICATES... --pop1 A --pop2 B --pop3 C --pop4 D [--alpha 0.05] [--relabel FILE] [options]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `REPLICATES` | text (paths) | required | Genotype prefixes or VCF files, one per simulated dataset. |
| `--pop1` | text | required | Population 1. |
| `--pop2` | text | required | Population 2. |
| `--pop3` | text | required | Population 3. |
| `--pop4` | text | required | Population 4. |
| `--alpha` | number (0 to 1) | 0.05 | A replicate rejects at p at or below it. |
| `--format` | `text`, `tsv`, `json` | `text` | `text`: the summary; `tsv`, `json`: one row per replicate. |
| `--output`, `-o` | text (path) | none | Also write one row per replicate as TSV here. |
| `--block-size` | number | 0.05 | jackknife block length in Morgans. |
| `--maxmiss` | number (0 to 1) | 0 | Keep SNPs missing in at most this fraction of the populations. |
| `--all-chromosomes` | flag | off | Use every chromosome, not only autosomes 1-22. |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | Which samples count one allele per SNP (see [SKILL.md](../SKILL.md#things-every-command-shares); `admixtools2` reproduces ADMIXTOOLS2-0007). |
| `--transversions-only` | flag | off | Use only transversions (drop A<->G and C<->T SNPs), for ancient DNA without UDG treatment, as admixtools2 `transitions = FALSE`. A SNP without two known bases is dropped. |
| `--transitions-only` | flag | off | Use only transitions, as admixtools2 `transversions = FALSE`. |
| `--minmaf` | number (0-0.5) | 0 | Drop SNPs whose minor allele frequency, the populations' alleles pooled, is below this (admixtools2 `minmaf`). |
| `--maxmaf` | number (0-0.5) | 0.5 | Drop SNPs whose pooled minor allele frequency is above this (admixtools2 `maxmaf`). |
| `--keep-snps` | text (path) | none | Use only the SNP IDs in this file (one per line, or a `.snp`/`.bim` file). As admixtools2 `keepsnps`, it overrides `--maxmiss`, the MAF bounds, `--transversions-only` and the autosome restriction (chromosome X SNPs listed are used); polymorphism and `--exclude-snps` still apply. |
| `--exclude-snps` | text (path) | none | Never use the SNP IDs in this file (popmix's own; admixtools2 has none). |
| `--keep-inds` | text (path) | none | Use only the sample IDs in this file (one per line, or a `.ind`/`.fam` file). |
| `--exclude-inds` | text (path) | none | Never use the sample IDs in this file (one per line, or a `.ind`/`.fam` file). |
| `--relabel` | text (path) | none | Give samples new populations: `ID POPULATION` per line, as admixtools2 `inds`/`pops`. Population names given to the command are those after relabelling. |
| `--relabel-pops` | text (path) | none | Rename populations: `OLD NEW` per line, applied after `--relabel` (merge two populations by giving them one name). |

## Inputs and outputs

Simulated datasets, which popmix only reads: make them with msprime (or any
simulator) in your own script and write each replicate as a VCF or an
EIGENSTRAT/PLINK dataset; msprime is GPL-3 and popmix never imports it.
Label VCF samples with `--relabel` (a file of `ID POPULATION` lines).
Writes a summary, or one row per replicate: `replicate`, `est`, `se`, `z`, `p` and `rejected`. In Python, `pm.calibrate_f4(replicates, pop1, pop2, pop3, pop4)`, with `attrs` as `calibrate_qpadm`'s.

```python
# not run: your own script, with msprime installed (popmix does not import it)
import msprime
for rep in range(100):
    demography = msprime.Demography()  # your model: populations, splits, admixture
    ...
    ts = msprime.sim_ancestry(samples={"A": 10, "B": 10, "C": 10, "D": 10}, demography=demography,
                              sequence_length=1e8, recombination_rate=1e-8, random_seed=rep + 1)
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=rep + 1)
    with open(f"sims/rep{rep}.vcf", "w") as f:
        ts.write_vcf(f)
# labels.txt: each VCF sample ID (tsk_0, tsk_1, ...) and its population
```

## Examples

```bash
# runs
# the example data given twice, to show the output: real replicates are simulations
popmix calibrate-f4 family family --pop1 French --pop2 Yoruba --pop3 Simulation --pop4 Simulation.1
```

```bash
# not run: needs simulated replicates
popmix calibrate-f4 sims/rep*.vcf --relabel labels.txt --pop1 A --pop2 B --pop3 C --pop4 D
```

## Reading the result

On a simulated tree ((A, B), (C, D)) f4 is 0, and the rate is its false-positive rate, which should be near `--alpha`. With simulated gene flow it is the power to detect it. The rate comes with a Wilson 95% interval: with 20 replicates a
true rate of 0.05 can show anywhere from 0 to 3 rejections, so simulate
enough replicates for the precision you need (100 gives about ±0.04).

## Caveats

- The rates hold for the simulated model, sample sizes, SNP density and
  missingness only; simulate what your data look like (ancient samples as
  pseudo-haploid and with missing calls, the SNP panel's ascertainment).
- f2 blocks are extracted for each replicate with the options given, as
  `--input` would.

## Reference

Simulation-based calibration of qpAdm (Harney et al. 2021) and of
f-statistics; msprime (Baumdicker et al. 2022) for the simulations.
