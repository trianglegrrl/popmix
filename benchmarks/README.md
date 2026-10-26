# Benchmarks against the reference tools

`run.py` runs popmix's command line and each analysis's reference tool on
the same input, several times, and records:

- wall time, CPU time and peak memory of every run (median, IQR,
  coefficient of variation; one unmeasured warm-up run first, so every
  measured run reads from the page cache);
- popmix at several thread counts, and whether its output is byte-identical
  at all of them;
- agreement of the two outputs at the tolerances of CLAUDE.md.

```bash
python benchmarks/run.py                       # every case, 5 repeats, 1/4/16 threads
python benchmarks/run.py --cases pca,kinship --repeats 3
python benchmarks/report.py benchmarks/results/<host>-<date>.json   # Markdown table
```

The input is AdmixTools 8 convertf's PACKEDANCESTRYMAP subset of AADR v66.p1
to 16 populations (473 individuals, 1,233,013 SNPs), made once in the
reference cache; the PCA and kinship cases use its 8 present-day
populations. Each case is in `cases.py`: the two command lines and the
comparison. The f2 cases run popmix with `--pseudohaploid admixtools2`, the
reference's ploidy rule (DEVIATIONS.md ADMIXTOOLS2-0007), and compare every
entry at atol 1e-15 + rtol 1e-10.

| Case | Reference |
|------|-----------|
| `convertf` | AdmixTools 8.0.2 convertf, packed to EIGENSTRAT |
| `extract_f2` | admixtools2 2.0.10 `f2_from_geno` |
| `qpadm` | admixtools2 2.0.10 `qpadm` on `f2_from_geno` blocks |
| `qpdstat` | admixtools2 2.0.10 `qpdstat` |
| `qpgraph` | admixtools2 2.0.10 `qpgraph` |
| `qpadm_rotate` | admixtools2 2.0.10 `qpadm_rotate` |
| `find_graphs` | admixtools2 2.0.10 `find_graphs` |
| `qpgraph_resample` | admixtools2 2.0.10 `qpgraph_resample_snps` |
| `classic_qpdstat`, `classic_qp3pop`, `classic_qpadm`, `classic_qpgraph` | ADMIXTOOLS 7.0.2 programs, compared at their printed precision (`classic_qpadm` differs by design: DEVIATIONS.md) |
| `pca` | EIGENSOFT 7.2.1 smartpca (full-precision build) |
| `kinship` | PLINK 2 `--make-king-table` |
| `read_kinship` | READv2 |
| `admixture` | ADMIXTURE 1.3.0 |
| `rolloff` | ADMIXTOOLS 7.0.2 rolloff (full-precision build) |
| `dates` | DATES v753 with its numbins fix |

Results are kept in `results/`, with the machine and the tool versions.

## Scaling

`scaling.py` measures how time and memory grow, popmix against the
reference where the reference finishes in reasonable time:

- populations: `extract_f2` on 4 to 497 AADR groups of 10 samples each
  (admixtools2 up to 256);
- SNPs: `extract_f2` on 1/16 to all of the 1,233,013 SNPs;
- samples: PCA (smartpca) and KING (PLINK 2) on 250 to 3,970 present-day
  samples;
- threads: 1 to 16 threads on `extract_f2` (128 populations), PCA and KING
  (2,000 samples);
- full AADR: `extract_f2` of the 16 populations from all 23,089 samples, from
  the packed copy and from the TGENO directly.

```bash
python benchmarks/scaling.py                   # every curve, 3 repeats
python benchmarks/scaling.py --curves snps,threads --repeats 1
```

Points already in the output file are skipped, so an interrupted run
resumes. `report_html.py` renders both result files as one page, with the
speed-ups, the agreement of every case and the scaling curves:

```bash
python benchmarks/report_html.py benchmarks/results/<host>-<date>.json \
    --scaling benchmarks/results/scaling-<host>-<date>.json -o report.html
```

## In CI

`smoke.py` times popmix's main commands on the small dataset shipped in the
repository (`tests/golden/dates/family_chr20_22`), with no reference tool;
the CI workflow (`.github/workflows/ci.yml`) writes its table to the job
summary, so a change in popmix's times shows from commit to commit:

```bash
python benchmarks/smoke.py --repeats 3
```
