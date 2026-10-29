# popmix

popmix does the population genetics analyses of admixtools2, ADMIXTOOLS,
EIGENSOFT smartpca, DATES, ROLLOFF, ADMIXTURE, PLINK 2's KING-robust and
READv2, in one Rust library with a Python package and a command line. It's
fast (often more than 10 times faster than the original tool), and every
analysis is tested against the original on real data.

## What it does

| Question | Analyses |
|----------|----------|
| How different are two populations? | f2, FST |
| Is a population admixed? | f3, qpAdm |
| Was there gene flow between populations? | f4, D-statistics, scans and matrices |
| What graph of splits and admixtures fits? | qpGraph, graph search |
| How many streams of ancestry are there? | qpWave |
| In what proportions? | qpAdm, ADMIXTURE |
| When did populations mix? | DATES, ROLLOFF |
| How are the samples structured and related? | PCA (as smartpca), KING-robust and READ kinship, site-frequency spectra |

It reads and writes VCF (read only), PLINK, EIGENSTRAT, PACKEDANCESTRYMAP,
TGENO (the AADR's format from v66 on), PED and ANCESTRYMAP, and merges
datasets as AdmixTools `mergeit` does.

## Install

```bash
pip install popmix          # library and command line
pip install popmix[plot]    # with plots (matplotlib)
conda install -c bioconda popmix
```

Wheels are built for Linux (x86_64 and ARM), macOS (Apple Silicon and
Intel) and Windows (x64), for Python 3.9 and later.

## A quick look

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

`aadr_v66p1` here is 27 populations of the AADR v66.p1, the set the tests
use (`tests/golden/aadr/build_aadr_v66p1.py` builds it). A typical analysis, from your own VCF to an admixture model and a
date:

```bash
popmix merge --file-a my_samples.vcf.gz \
  --file-b /references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB \
  --output merged --output-format plink
popmix qpadm -i merged --target MySample \
  --left French,Han --right Yoruba,Mbuti,Papuan
popmix dates -i merged --admixed MySample --sources French,Han
```

And in Python:

```python
import popmix as pm

f2 = pm.extract_f2("merged", populations=["MySample", "French", "Han", "Yoruba", "Mbuti", "Papuan"])
fit = pm.qpadm(f2, target="MySample", left=["French", "Han"], right=["Yoruba", "Mbuti", "Papuan"])
print(fit.weights_dict(), fit.p)
print(fit.assess())  # whether the model fits, and why or why not
```

[docs/guide.md](docs/guide.md) walks through every command with real
output, the Python API, the file formats and the error messages. Every
command also documents its options under `popmix <command> --help`.

## How much to trust it

Each analysis is checked against its original tool on the AADR, at
tolerances tight enough that a wrong formula can't hide (f-statistics to a
relative 1e-10, standard errors and qpAdm weights to 1e-6):

| Analysis | Checked against | Agreement |
|----------|-----------------|-----------|
| f2, f3, f4, D, FST, scans | admixtools2 2.0.10 | 1e-10 (SEs 1e-6) |
| qpWave, qpAdm, resampling | admixtools2 2.0.10 | chi-squared, SEs and weights 1e-6, p 1e-4 |
| qpGraph and graph search | admixtools2 2.0.10 | scores 3e-12; same graphs and optima |
| PCA | EIGENSOFT 7.2.1 smartpca | eigenvalues 1e-9 |
| DATES | DATES v753 | printed precision |
| ROLLOFF | ADMIXTOOLS 7.0.2 | curves 1e-9 |
| Kinship | PLINK 2 (KING-robust), READv2 | 1e-10 |
| ADMIXTURE | ADMIXTURE 1.3.0 | log-likelihood at least ADMIXTURE's |
| SFS | admixtools2, dadi 2.4.4 | exact, 1e-10 |

[TESTING.md](TESTING.md) lists every check. Where an original tool has a
bug, popmix fixes it and a `--compat` option reproduces the original;
[DEVIATIONS.md](DEVIATIONS.md) describes each one, with the runs that show
it.

The statistics are only as good as the model behind them. A qpAdm model
that isn't rejected isn't proven (popmix lists the simpler models that fit
too), qpWave's rank is a lower bound, and one DATES date doesn't mean one
pulse of admixture. The guide's [caveats](docs/guide.md#caveats) cover
these and a few more.

## Speed

On the same 473 AADR v66.p1 individuals (1.23M SNPs), median of 5 runs on
a 16-thread laptop:

| Analysis | Original | popmix |
|----------|----------|--------|
| f2 blocks, 16 populations | admixtools2: 36.9 s, 7.9 GB | 0.9 s, 0.7 GB |
| qpAdm from genotypes | admixtools2: 21.5 s | 1.0 s |
| PCA, 10 PCs | smartpca: 83.2 s | 4.3 s |
| ROLLOFF | ADMIXTOOLS 7.0.2: 4.3 min | 15.7 s |
| DATES | DATES v753: 30.7 s | 2.3 s |
| READ kinship | READv2: 28.3 s | 0.8 s |
| KING-robust kinship | PLINK 2: 1.2 s | 1.1 s |

popmix is faster than the original in every benchmark, though only just
for KING (at 3,970 samples PLINK 2 is 7% faster). The output is identical at 1,
4 and 16 threads. [benchmarks/](benchmarks/) has the full results, with
every analysis, memory and scaling to 23,089 samples.

## Agent skill

popmix ships a skill for Claude Code and other agents
([skills/popmix/](skills/popmix/)): which command answers which question,
and every option with its default, examples and caveats. Its examples run
in the test suite, so it can't drift from the command line.

```bash
popmix install-skill          # into ~/.claude/skills/popmix
popmix install-skill --force  # after upgrading popmix
```

## Building from source

```bash
git clone https://github.com/trianglegrrl/popmix
cd popmix
python -m venv .venv && source .venv/bin/activate
pip install maturin pytest hypothesis pyyaml numpy pandas click
maturin develop --release -m crates/popmix-python/Cargo.toml
cargo test --workspace && pytest tests/python
```

## API stability

The Python package (the names in `popmix.__all__`, their parameters and
results) and the command line (commands, options and output formats) follow
semantic versioning from 1.0: a 1.x release adds to them but doesn't break
them. The Rust crates (`popmix-core`, `popmix-io`, `popmix-python`) are
popmix's internals. You can build on them, but their Rust API has no
stability promise, and it may change in any release, 1.x included. An analysis may still correct a result against
its original tool in a 1.x release; DEVIATIONS.md lists each such change.

## Citation

Please cite popmix as [CITATION.cff](CITATION.cff) gives it (GitHub's "Cite
this repository" shows it as APA and BibTeX), and the papers of the methods
you use. popmix reimplements other people's methods, and they deserve the
credit:

| Analysis | Papers |
|----------|--------|
| f2, f3, f4, D, FST | Reich et al. 2009; Green et al. 2010; Patterson et al. 2012; Maier et al. 2023 |
| qpWave, qpAdm | Reich et al. 2012; Haak et al. 2015; Harney et al. 2021 |
| Admixture graphs (qpGraph) | Patterson et al. 2012; Maier et al. 2023 |
| Block jackknife | Busing et al. 1999 |
| PCA, Tracy-Widom | Patterson, Price and Reich 2006; Price et al. 2006; Tracy and Widom 1994 |
| ROLLOFF | Moorjani et al. 2011 |
| DATES | Narasimhan et al. 2019; Chintalapati et al. 2022 |
| ADMIXTURE | Alexander et al. 2009; Varadhan and Roland 2008 (SQUAREM) |
| KING-robust kinship | Manichaikul et al. 2010 |
| READ kinship | Monroy Kuhn et al. 2018; Alaçamlı et al. 2024 |
| Site-frequency spectra | Gutenkunst et al. 2009 |
| qpAdm and f4 calibration | Baumdicker et al. 2022 (msprime) |
| File formats, conversion, merge | Chang et al. 2015; Danecek et al. 2021 |

The admixture-graph fits match R's `optim` L-BFGS-B and its random starts
R's Mersenne Twister (Byrd et al. 1995; Zhu et al. 1997; Morales and
Nocedal 2011; Moré and Thuente 1994; Matsumoto and Nishimura 1998). The
tests run on the AADR (Mallick et al. 2024). `CITATION.cff` lists all of
these with their DOIs, and every original tool with the version the tests
ran; [PROVENANCE.md](PROVENANCE.md) records what each module was written
from and how popmix was built.

## License

MIT (LICENSE). The parts derived from ADMIXTOOLS, EIGENSOFT, L-BFGS-B and
the Mersenne Twister reference code keep their BSD 3-Clause licences
(NOTICE), and THIRDPARTY.yml holds the licences of the Rust crates compiled
into popmix.

popmix is independent: it is not affiliated with or endorsed by the authors
of admixtools2, ADMIXTOOLS, EIGENSOFT, DATES, ROLLOFF, ADMIXTURE, PLINK,
KING or READv2. Their names identify the tools whose results it matches.
