# The family data and the DATES references

## The family data

`family_chr20_22.{geno,snp,ind}` (in the repository) and
`example_data/family_packed.{geno,snp,ind}` (every autosome, 42 MB, built
locally and gitignored) are PACKEDANCESTRYMAP files of the same 50 samples:

| Population | Samples | What |
|------------|---------|------|
| `French` | 24 | AADR v66.p1 `French.DG`: HGDP00511.DG to HGDP00535.DG, but HGDP00532 |
| `Yoruba` | 21 | AADR v66.p1 `Yoruba.DG`: HGDP00920.DG, HGDP00924.DG to HGDP00926.DG, HGDP00928.DG to HGDP00944.DG |
| `Simulation` | 1 | NA00000, simulated: 20% French, 80% Yoruba, 50 generations ago |
| `Simulation.1` | 4 | NA00001 to NA00004, simulated as NA00000 |

`family_truth.json` lists every sample ID. The French and Yoruba genomes
are AADR v66.p1's own (1240K public release, doi:10.7910/DVN/FFIDCW,
licensed CC0), diploid calls from whole-genome sequence (the `.DG` data
type), relabelled `French` and `Yoruba`; sex is AADR's, `U` for the
simulated samples.

SNPs: every third autosomal SNP of the 1240K panel (rows 0, 3, 6, ... of
its `.snp`): 383,547 in `family_packed`, of which the 21,175 on chromosomes
20 to 22 are `family_chr20_22`. Positions, genetic map and alleles are
AADR's.

### The simulation

Each simulated sample is two haplotypes, made chromosome by chromosome:

- ancestry changes along a haplotype at the crossovers since the admixture,
  a Poisson process of 50 per Morgan on the AADR genetic map (the `.snp`
  file's genetic position), from the chromosome's first SNP to its last;
- each segment between crossovers is French with probability 0.2, else
  Yoruba, independently, so ancestry LD decays as exp(-50 d), the model
  DATES and ROLLOFF fit;
- a segment's alleles come from one source genome of its ancestry, chosen
  at random, never a reference sample: for French, HGDP00536.DG to
  HGDP00539.DG (AADR `French.DG`); for Yoruba, 30 genomes of AADR's
  `YRI.DG` (the 1000 Genomes Yoruba in Ibadan, and HGDP00927.DG, which
  AADR labels YRI), the first by sample ID that are unrelated to every HGDP
  Yoruba and to each other (KING-robust kinship below 2^-4.5, third
  degree). `family_truth.json` lists them. The sources are unphased, so at
  each SNP the haplotype takes one of the source's two alleles at random; a
  missing call stays missing;
- the genotype is the sum of the two haplotypes, missing if either is.

The Yoruba pool is large on purpose. Every simulated sample draws from the
same sources, so the sources' own sampling noise acts as drift shared by the
whole simulated population: with 4 French and 6 Yoruba source genomes,
f3(Simulation.1; French, Yoruba) came out at zero, not negative as for an
admixed population; with 30 Yoruba sources it is clearly negative
(`tests/python/test_family_data.py`).

The other French.DG and Yoruba.DG genomes of v66.p1 (B_French-3, S_French-1,
S_French-2, B_Yoruba-3, S_Yoruba-1, S_Yoruba-2: SGDP and Meyer et al. 2012)
are HGDP individuals sequenced again (their calls agree with one HGDP
genome's at more than 99.9% of SNPs) and are left out. No two of the
genomes used, references or sources, are related (KING-robust kinship at
most 0.013), and none of the kinship goldens and tests needs related pairs,
so none are simulated.

numpy's `default_rng(20261029)` draws everything in a fixed order.

### The truth

`family_truth.json` records the parameters, the sources and, for each
haplotype, its tracts (about 1,800 genome-wide) and the French share of its
genetic length. The realized French ancestry (mean of the two haplotypes):

| Sample | Genome | Chromosomes 20-22 |
|--------|--------|-------------------|
| NA00000 (Simulation) | 0.191 | 0.146 |
| NA00001 | 0.197 | 0.173 |
| NA00002 | 0.227 | 0.228 |
| NA00003 | 0.202 | 0.230 |
| NA00004 | 0.208 | 0.179 |


### Scripts

- `build_family.py` builds both datasets and `family_truth.json` from the
  AADR v66.p1 subset that `tests/golden/aadr/build_aadr_v66p1.py` writes in
  the reference cache (checked against its manifest), with AdmixTools 8
  convertf writing the PACKEDANCESTRYMAP files.
  `tests/golden/test_family_build.py` (`--run-slow`) reruns it and checks
  the output is byte-identical; `tests/python/test_family_data.py` checks
  the composition and that supervised ADMIXTURE finds the simulated
  ancestry.
- `build_dates.sh` builds DATES v753 (with `dates_numbins_fix.patch`);
  `generate_dates_reference.sh` runs it on these data and writes
  `v753/<config>/`: the parameter file, the LD table (`dates.out`), the
  full-data fit (`dates.flog`, `dates.fit`), the per-chromosome jackknife
  (`dates.jin`) and its result (`dates.jout`). `chr20_22_*` read
  `family_chr20_22`, `simulation_*` read `example_data/family_packed`;
  `_q0` is DATES's exact mode, `_q10` its FFT mesh.
- `generate_yamnaya_subset.sh` writes `yamnaya_chr21_22`, AADR v66.p1
  samples for a source of mixed ploidy (DATES-0004).

The full-genome data are not in the repository: run `build_family.py`
before the slow DATES goldens (`simulation_q0`, `simulation_q10`), which
skip without them.
