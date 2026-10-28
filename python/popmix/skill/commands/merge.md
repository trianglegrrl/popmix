# popmix merge

Merge two genotype datasets into one: every sample of both, on the SNPs they
share, with alleles reconciled as AdmixTools mergeit does. Typical use: add
your own samples (VCF or PLINK) to a reference panel such as the AADR.

## Synopsis

```text
popmix merge --file-a PREFIX --file-b PREFIX --output PREFIX --output-format FORMAT
             [--input-format-a FORMAT] [--input-format-b FORMAT]
```

## Arguments and options

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--file-a`, `-a` | text (path) | required | First dataset (prefix, or full VCF path). Its samples come first, and file B's genotypes are recoded to its alleles. |
| `--file-b`, `-b` | text (path) | required | Second dataset (prefix, or full VCF path). |
| `--output`, `-o` | text (path) | required | Output prefix. |
| `--output-format`, `-of` | `ancestrymap`, `eigenstrat`, `ped`, `packedped`, `plink`, `packedancestrymap`, `tgeno` | required | Output format (no VCF). |
| `--input-format-a`, `-ifa` | `ancestrymap`, `eigenstrat`, `ped`, `packedped`, `plink`, `packedancestrymap`, `tgeno`, `vcf` | none (auto-detected) | Format of file A. |
| `--input-format-b`, `-ifb` | `ancestrymap`, `eigenstrat`, `ped`, `packedped`, `plink`, `packedancestrymap`, `tgeno`, `vcf` | none (auto-detected) | Format of file B. |
| `--strandcheck`/`--no-strandcheck` | flag | `--strandcheck` | Reconcile strands: complemented alleles are the same SNP on the other strand, and A/T and C/G SNPs are left out because their strand cannot be told. `--no-strandcheck` when both files are on one strand (two releases of one panel, such as two AADR 1240K releases): A/T and C/G SNPs are kept, and complemented alleles count as a mismatch, as mergeit `strandcheck: NO`. |

## Inputs and outputs

Writes the merged dataset under `--output` and prints a report:

```text
Merge complete:
  SNPs in file A: 5
  SNPs in file B: 4
  SNPs in common: 3
  Samples in file A: 2
  Samples in file B: 1
  Total samples: 3
```

with, when they occur, `Alleles swapped (recoded)`, `On the other strand
(complemented)` and one `Left out, <reason>: N` line per reason.

## Examples

```bash
# runs
popmix merge -a merge_a -b merge_b -o merged -of eigenstrat
popmix populations merged --format tsv
popmix merge -a merge_a -b merge_b -o merged_plink -of plink
```

```bash
# not run: needs the AADR and your own VCF
popmix merge -a aadr_v66p1 -b my_samples.vcf.gz \
  -o panel_plus_mine -of packedancestrymap
```

## Reading the result

When the two datasets share less than 90% of the smaller one's SNPs they
were ascertained on different panels (the AADR's 1240K and Human Origins,
say), and popmix warns: the merged data keep only the shared SNPs, a panel
of their own, so compare f-statistics only among populations taken from
the merged data, never with f2 blocks extracted from either input.

`SNPs in common` is the number of SNPs written. SNPs are paired by
chromosome and position. For each pair: alleles in the other order are
recoded; alleles that are the complement (the other strand) are flipped; the
missing allele of a SNP monomorphic in one file is taken from the other.
Left out and counted: A/T and C/G SNPs (strand unknown), allele mismatches
(a third allele), and positions where file B has several SNPs and none
with file A's name.

## Caveats

- SNPs are paired by position, not by name, unlike mergeit and PLINK, so
  the same SNP under an rs ID and a chr:pos name is merged; both files must
  use the same genome build
  ([mergeit interface differences](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#mergeit-admixtools-802)).
- Dropping A/T and C/G SNPs is mergeit's `strandcheck`.
- A sample ID present in both files is an error (`Duplicate sample IDs found
  in both files`); rename or remove it in one file first.
- Each sample keeps the sex its input file gives.
- Genetic positions are written with six decimals.

## Reference

AdmixTools `mergeit` (8.0.2), checked against it in
tests/golden/test_merge_mergeit.py.
