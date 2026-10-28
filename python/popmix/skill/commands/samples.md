# popmix samples

Each sample's coverage: SNPs called, the share of calls that are
heterozygous, and the ploidy each rule finds. Use it before an analysis
with ancient or low-coverage samples, to decide which to project in a PCA,
whether KING-robust or READ suits the kinship, and how much a ploidy call
rests on.

## Synopsis

```text
popmix samples --input PREFIX [--populations P1,P2,...] [--all-chromosomes] [--format tsv|json] [-o FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix, or the full path of a VCF. |
| `--populations` | text (comma-separated) | none (all but `Ignore`) | Populations whose samples are listed. |
| `--all-chromosomes` | flag | off | Count every chromosome, not only autosomes 1-22. |
| `--format` | `tsv`, `json` | `tsv` | Output format. |
| `--output`, `-o` | text (path) | none (stdout) | Write the table here. |

## Inputs and outputs

One row per sample, in file order: `id`, `pop`, `snps_called`,
`called_fraction` (of the SNPs read), `heterozygosity` (the share of calls
that are heterozygous), `ploidy` (1 pseudo-haploid, 2 diploid, as popmix
finds it: no heterozygous call at any SNP) and `ploidy_admixtools2` (no
heterozygous call among the first 1,000 SNPs, as admixtools2). In Python,
`pm.sample_coverage(prefix, populations)`.

## Examples

```bash
# runs
popmix samples -i family
popmix samples -i family6 --populations French_1,Simulation.1 --format json
```

## Reading the result

- A pseudo-haploid sample (one allele read per SNP, coded 0 or 2) has
  heterozygosity 0 and ploidy 1; f-statistics count one allele for it.
  Where the two ploidy columns differ, admixtools2's rule has seen too few
  SNPs (ADMIXTOOLS2-0007).
- Samples with a low `called_fraction` are best projected onto a PCA built
  from well-covered ones (`pca --project-pops ... --lsqproject`); `pca`
  warns when a PCA sample misses more than half the SNPs.
- KING-robust needs heterozygous calls and many shared SNPs (it warns
  under 10,000); for pseudo-haploid or low-coverage samples use
  `kinship --method read`.

## Caveats

- Coverage here is calls in the genotype file, not sequencing depth: a
  pseudo-haploid call at 0.1x and a diploid call at 30x both count once.
- Heterozygosity depends on the SNP panel's ascertainment; compare it
  between samples on one panel, not with published genome-wide values.

## Reference

None; counts of the genotype file's calls.
