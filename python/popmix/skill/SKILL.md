---
name: popmix
description: Use when the user wants population-genetics analyses of genotype data with the popmix command line - f2/f3/f4 and D-statistics, qpAdm/qpWave admixture models, qpGraph admixture graphs, admixture dating (DATES, ROLLOFF), smartpca PCA, ADMIXTURE ancestry proportions, KING/READ kinship, site-frequency spectra, or converting and merging EIGENSTRAT, PACKEDANCESTRYMAP, TGENO, PLINK and VCF files - or asks how to do in popmix what admixtools2, ADMIXTOOLS, smartpca, DATES, ADMIXTURE, PLINK 2 KING or READv2 do.
---

# popmix

This skill documents popmix 1.0.0.

popmix is an MIT-licensed reimplementation, in Rust with a Python package and
a `popmix` command line, of the standard tools for ancient and present-day
population genetics: admixtools2 (f-statistics, qpAdm, qpWave, qpGraph,
graph search), ADMIXTOOLS (ROLLOFF, convertf, mergeit), DATES, EIGENSOFT
smartpca, ADMIXTURE, PLINK 2's KING-robust kinship, READv2 and admixtools2 /
dadi site-frequency spectra. Every analysis is checked against its reference
tool on real data (AADR v66.p1); where the reference has a bug popmix
corrects it by default and reproduces it with a `--compat` flag (see
[Deviations](#deviations-and-compat-flags)).

Check what is installed with `popmix --version`; if it is not this version,
some options here may differ, and `popmix <command> --help` is authoritative.
Reinstall the matching skill with `popmix install-skill --force`.

## Which command for which question

| Question | Command | Reference tool |
|----------|---------|----------------|
| What is in this dataset? Which populations, how many samples? | [info](commands/info.md), [populations](commands/populations.md) | |
| Convert between formats (VCF to PLINK, PLINK to EIGENSTRAT, ...) | [convertf](commands/convertf.md) | AdmixTools convertf |
| Combine my samples with a reference panel | [merge](commands/merge.md) | AdmixTools mergeit |
| Precompute f2 blocks once for many analyses | [extract-f2](commands/extract-f2.md) | admixtools2 `extract_f2()` |
| How much drift separates A and B? | [f2](commands/f2.md) | admixtools2 `f2()` |
| Is A admixed between B and C? (negative f3) | [f3](commands/f3.md), [qp3pop](commands/qp3pop.md) | admixtools2 `f3()`, `qp3pop()` |
| Gene flow between branches (ABBA-BABA, f4, D) | [f4](commands/f4.md), [qpdstat](commands/qpdstat.md) | admixtools2 `f4()`, `qpdstat()` |
| Which populations share most drift with X? | [outgroup-f3](commands/outgroup-f3.md) | admixtools2 `qp3pop(outgroupmode = TRUE)` |
| D or f4 for every pair of a population list, as a heatmap | [dstat-matrix](commands/dstat-matrix.md) | admixtools2 `qpdstat()` |
| Is target T a mixture of sources S1..Sn, and in which proportions? | [qpadm](commands/qpadm.md) | admixtools2 `qpadm()` |
| How many ancestry streams relate these populations? | [qpwave](commands/qpwave.md) | admixtools2 `qpwave()` |
| Does this admixture graph fit? | [qpgraph](commands/qpgraph.md) | admixtools2 `qpgraph()` |
| Which admixture graph fits best? | [find-graphs](commands/find-graphs.md) | admixtools2 `find_graphs()`, `generate_all_graphs()` |
| When did the admixture happen? | [dates](commands/dates.md) (recommended), [rolloff](commands/rolloff.md) | DATES v753, ADMIXTOOLS rolloff |
| PCA, with ancient samples projected | [pca](commands/pca.md) | EIGENSOFT smartpca |
| Ancestry proportions with K clusters | [admixture](commands/admixture.md) | ADMIXTURE 1.3.0 |
| Are these samples related? | [kinship](commands/kinship.md) | PLINK 2 `--make-king-table`, READv2 |
| Site-frequency spectrum, Tajima's D | [sfs](commands/sfs.md) | dadi, admixtools2 `joint_sfs()` |
| Install this skill for Claude Code | [install-skill](commands/install-skill.md) | |

A typical f-statistics session: `populations` to find the exact names, then
either pass `--input PREFIX` to each command (f2 blocks are extracted for
that command's populations), or run `extract-f2` once for every population
of the project and pass `--f2-dir` to the rest, which is faster and puts
every statistic on the same SNP set, as admixtools2 does with precomputed f2.

## Command index

Each file gives the synopsis, every option with type, default and meaning,
inputs and outputs, examples, how to read the result, caveats and the
reference tool.

- [admixture](commands/admixture.md): ancestry proportions with K ancestral populations
- [convertf](commands/convertf.md): convert between genotype formats
- [dates](commands/dates.md): admixture date from weighted LD (DATES)
- [dstat-matrix](commands/dstat-matrix.md): D or f4 for every pair of populations
- [extract-f2](commands/extract-f2.md): compute and save f2 blocks
- [f2](commands/f2.md): f2 between two populations
- [f3](commands/f3.md): one f3 admixture test
- [f4](commands/f4.md): one f4 statistic
- [find-graphs](commands/find-graphs.md): search for the best admixture graphs
- [info](commands/info.md): dataset summary
- [install-skill](commands/install-skill.md): install this skill
- [kinship](commands/kinship.md): KING-robust or READ-style kinship of sample pairs
- [merge](commands/merge.md): merge two datasets
- [outgroup-f3](commands/outgroup-f3.md): outgroup-f3 matrix
- [pca](commands/pca.md): smartpca PCA with projection
- [populations](commands/populations.md): populations and sample counts
- [qp3pop](commands/qp3pop.md): many f3 statistics
- [qpadm](commands/qpadm.md): admixture weights of a target
- [qpadm-rotate](commands/qpadm-rotate.md): every qpAdm model of a rotation, one row each
- [qpdstat](commands/qpdstat.md): many f4 or D statistics
- [qpgraph](commands/qpgraph.md): fit an admixture graph
- [qpwave](commands/qpwave.md): rank test of left against right populations
- [rolloff](commands/rolloff.md): admixture date by ROLLOFF
- [sfs](commands/sfs.md): site-frequency spectra

## Global options

Global options go before the command: `popmix --threads 8 qpadm ...`.

| Option | Type | Default | Meaning |
|--------|------|---------|---------|
| `--threads`, `-t` | integer | 4 | Worker threads (rayon pool) for every command; at least 1. Results are byte-identical at any thread count. |
| `--debug` | flag | off | On an error, show the Python traceback instead of the one-line `Error: ...` message. |
| `--version` | flag | off | Print `popmix, version X` and exit. |

`popmix --help` lists the commands and `popmix <command> --help` its options.

## Input formats

Most commands take a genotype **prefix**: the path without the extension,
shared by the dataset's files. A VCF is given by its full path, with the
extension. The format is detected from the files:

| Format | Files | Notes |
|--------|-------|-------|
| PACKEDANCESTRYMAP | `.geno` (binary, `GENO` header), `.snp`, `.ind` | Packed EIGENSTRAT, one record per SNP. Streamed SNP by SNP by `extract-f2` and the f-statistic commands, so it is the format to use for large data. |
| TGENO | `.geno` (binary, `TGENO` header), `.snp`, `.ind` | Packed, one record per individual (AADR v66 and later). |
| EIGENSTRAT | `.geno` (text, one digit per sample: 0, 1, 2 or 9 for missing), `.snp`, `.ind` | `.snp`: id, chromosome, genetic position (Morgans), physical position, two alleles. `.ind`: sample, sex, population. |
| PLINK (PACKEDPED) | `.bed`, `.bim`, `.fam` | The population is the family ID (first `.fam` column). Genetic positions are read in Morgans; a `.bim` whose largest position is over 10 is taken as centimorgans. |
| PED | `.ped`, `.map`, `.pedind` | Text PLINK. |
| ANCESTRYMAP | `.ancestrymapgeno`, `.snp`, `.ind` | Sparse text. |
| VCF | `.vcf` or `.vcf.gz` (full path) | Input only. Biallelic SNPs only (others are skipped and counted); FILTER is not read; each sample is its own population, so merge with a labeled panel or convert and edit the `.ind`. |

Population labels come from the `.ind` (or `.fam`); relabeling samples, or
setting a sample's population to `Ignore`, is done by editing that file.
`pca`, `kinship`, `admixture` and `extract-f2` leave out `Ignore` when they
default to all populations.

The genotype counts the first allele of the `.snp` (EIGENSTRAT) or A1
(PLINK). Analyses that use genetic distance (f2 blocks, DATES, ROLLOFF)
need the genetic map column; with a map of all zeros, f2 blocks fall back to
2 Mb of physical distance, as admixtools2.

## Things every command shares

- **Population names** are checked before any work: an unknown name stops
  with "Population 'X' not found", suggestions and the command to list them.
  Population lists are comma-separated with no spaces needed
  (`--left French,Yoruba`); quote names with shell characters.
- **Errors**: bad input gives one line `Error: ...` on stderr and exit code
  1 (Click usage errors exit 2); `--debug` shows the traceback.
- **f2 sources**: the f-statistic commands take `--input PREFIX` (genotypes)
  or `--f2-dir DIR` (from `extract-f2`), never both.
- **Pseudo-haploid samples** (ancient samples with one allele called per
  SNP, coded 0 or 2) count one allele per SNP in f-statistics, f2 blocks,
  SFS and DATES. Which samples are pseudo-haploid is set by
  `--pseudohaploid` where a command has it: `detect` (default) takes those
  with no heterozygous call at any SNP; `admixtools2` takes those with none
  among the first 1,000 SNPs, which also catches some diploid genomes
  ([ADMIXTOOLS2-0007](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md#admixtools2-0007-diploid-genomes-taken-as-pseudo-haploid)),
  to reproduce admixtools2 exactly; `diploid` counts every sample as
  diploid. PCA, ADMIXTURE and KING read genotypes as called, as their
  reference tools do.
- **Output**: single statistics print text (or `--format json` / `tsv`);
  scans print a TSV table (or JSON) to stdout or `-o FILE`; PCA, ADMIXTURE
  and ROLLOFF write files under `--output-prefix`. Progress messages of
  `dates` and `kinship` go to stderr, so stdout stays parseable.
- **Provenance**: a command that writes files leaves
  `OUTPUT.provenance.json` beside its main output (`PREFIX.provenance.json`
  for an output prefix, `DIR/provenance.json` in an f2 directory): popmix
  version, command line, every parameter with its default filled in, and
  the size and modification time of each input. Output on stdout gets none.
- **Plots** (`--plot`, `--heatmap`, `--plot-fit`) need matplotlib
  (`pip install popmix[plot]`); the format follows the file extension.

## Deviations and compat flags

[DEVIATIONS.md](https://github.com/trianglegrrl/popmix/blob/main/DEVIATIONS.md)
(in the popmix repository) lists every place where popmix
corrects its reference tool. The corrections are the default; these
switches reproduce the reference exactly, for checking against it:

| Switch | Command | Reproduces |
|--------|---------|------------|
| `--compat dates-v753` | [dates](commands/dates.md) | DATES-0001 to DATES-0004 |
| `--compat eigensoft-7.2.1` | [pca](commands/pca.md) | EIGENSOFT-0001 (Tracy-Widom p-values beyond the table) |
| `--compat readv2-2.01` | [kinship](commands/kinship.md) `--method read` | READV2-0001 to READV2-0003 |
| `--compat admixtools2-2.0.10` | [sfs](commands/sfs.md) `--joint` | ADMIXTOOLS2-0003 |
| `--pseudohaploid admixtools2` | f-statistic commands, `extract-f2` | ADMIXTOOLS2-0007 |

## Example data

Every example marked `# runs` in these files is executed by popmix's test
suite (tests/python/test_skill.py) in a directory holding copies of the
small data bundled in the repository, so the examples cannot go stale:

- `family.geno/.snp/.ind`: tests/golden/dates/family_chr20_22, PACKEDANCESTRYMAP,
  54 samples and 21,162 SNPs on chromosomes 20-22: French (28), Yoruba
  (21) and simulated admixed samples, Simulation (1) and Simulation.1 (4),
  20% French and 80% Yoruba ancestry, 50 generations ago.
- `family6.*`: the same genotypes with French and Yoruba each split into
  three groups by sample order (French_1..3, Yoruba_1..3), so qpAdm and
  qpWave have enough populations. A toy labeling, for the examples only.
- From tests/fixtures: `small_eigenstrat.*` and `small_tgeno.*` (3 samples,
  5 SNPs, populations PopA and PopB), `small_vcf.vcf`, and `merge_a.*`,
  `merge_b.*` (two EIGENSTRAT sets that share 3 SNPs).

Examples marked `# not run: ...` need large data such as the AADR
(`aadr_v66p1` in these files: AADR v66.p1 with each population labelled by
its data type, as `French.DG`); adapt the path, and check the population
names, which are illustrative, with `popmix populations PREFIX --filter`.
Use AADR v66.p1 or later, never v62.0 or v66.0: those hold the 161 genomes
of Jacobs et al. 2019 that v66.p1 withdrew.

```bash
# runs
popmix --version
popmix info family
popmix populations family --format tsv
```
