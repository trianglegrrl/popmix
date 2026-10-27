# popmix sfs

Site-frequency spectrum of one or more populations: projected to chosen
sample sizes as dadi does (the right spectrum with missing data), or the
joint spectrum of admixtools2 `joint_sfs()`. For one population, also the
number of segregating sites, pi, Watterson's theta and Tajima's D.

## Synopsis

```text
popmix sfs --input PREFIX --populations A[,B,...] [--projection N1,N2,...] [--folded]
           [--joint [--compat admixtools2-2.0.10]] [--all-chromosomes]
           [--pseudohaploid RULE] [--format tsv|json] [-o FILE]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--input`, `-i` | text (path) | required | Genotype prefix. |
| `--populations` | text (comma-separated) | required | Populations, one axis of the spectrum each. |
| `--projection` | text (comma-separated integers) | none (full sample sizes) | Allele counts to project each population down to (one per population), by hypergeometric sampling, as dadi `Spectrum.from_data_dict`. SNPs with fewer observed alleles are left out; without a projection, SNPs with any missing data are left out. Diploid samples carry 2 alleles, pseudo-haploid samples 1. |
| `--folded` | flag | off | Fold the spectrum (minor-allele counts), as dadi `fold()`. Use it unless the counted allele is the derived one. |
| `--joint` | flag | off | admixtools2 `joint_sfs()` instead: counts among the observed alleles, no projection, every chromosome; `--projection` and `--folded` are refused, and `--all-chromosomes` is implied. |
| `--compat` | `admixtools2-2.0.10` | none (corrected) | Only with `--joint`: reproduce admixtools2 exactly, which drops SNPs whose count is not an exact integer in floating point (ADMIXTOOLS2-0003). |
| `--all-chromosomes` | flag | off | Use every SNP, not only autosomes (`--joint` always uses every SNP, as admixtools2). |
| `--pseudohaploid` | `detect`, `admixtools2`, `diploid` | `detect` | Which samples count one allele per SNP (see [SKILL.md](../SKILL.md#things-every-command-shares); `admixtools2` reproduces ADMIXTOOLS2-0007). |
| `--format` | `tsv`, `json` | `tsv` | Table format: tab-separated with a header, or JSON records. |
| `--output`, `-o` | text (path) | none (stdout) | Write the spectrum table here instead of stdout. |

## Inputs and outputs

A TSV table, one row per spectrum entry: one column per population (the
allele count in it) and `n`, the number of SNPs (expected, after
projection); the projected spectrum also has `masked` (dadi's mask: the
corners, allele absent or fixed everywhere, and with `--folded` the entries
folded away). For one population, unfolded and not `--joint`, four lines
go to stderr:

```text
S	17667.4
pi	6276.47
watterson_theta	4979.9
tajima_d	1.09285
```

## Examples

```bash
# runs
popmix sfs -i family --populations Yoruba --projection 20 -o yoruba_sfs.tsv
popmix sfs -i family --populations French,Yoruba --projection 10,10 --folded -o joint_folded.tsv
popmix sfs -i family --populations French,Yoruba --joint -o joint.tsv
popmix sfs -i family --populations French,Yoruba --joint --compat admixtools2-2.0.10 -o joint_at2.tsv
```

```bash
# not run: needs the AADR
popmix sfs -i aadr_v66p1 --populations Yoruba.DG --projection 20 --folded
```

## Reading the result

`n` at count i is the number of SNPs (expected, projected) with i copies of
the counted allele; ignore masked rows. Tajima's D near 0 fits a neutral
constant-size population; positive values point to a contraction or
balancing selection, negative to expansion. On SNP-array data (such as the
1240k panel), ascertainment distorts the spectrum and these statistics:
rare variants are missing, so D is biased upward.

## Caveats

- The counted allele is the genotype file's (EIGENSTRAT first allele, PLINK
  A1), so the unfolded spectrum is polarized only if that allele is the
  derived one; otherwise use `--folded`.
- `--joint` counts a SNP with missing data at its observed count, as
  admixtools2; prefer a projection for data with missing calls.
- Heterozygous calls in pseudo-haploid samples (half alleles) are left out.
- X and Y are counted with autosomal ploidy under `--joint`, as admixtools2.

## Reference

dadi 2.4.4 (`Spectrum.from_data_dict`, `fold`, `Tajima_D`) to 1e-10, and
admixtools2 `joint_sfs()` (2.0.10) exactly (with `--compat`).
