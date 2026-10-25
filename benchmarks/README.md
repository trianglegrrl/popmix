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
comparison.

| Case | Reference |
|------|-----------|
| `convertf` | AdmixTools 8.0.2 convertf, packed to EIGENSTRAT |
| `extract_f2` | admixtools2 2.0.10 `f2_from_geno` |
| `qpadm` | admixtools2 2.0.10 `qpadm` on `f2_from_geno` blocks |
| `pca` | EIGENSOFT 7.2.1 smartpca (full-precision build) |
| `kinship` | PLINK 2 `--make-king-table` |

Results are kept in `results/`, with the machine and the tool versions.
