# Popmix nf-core Modules

Nextflow DSL2 modules for running popmix population genetics analyses in nf-core pipelines.

## Available Modules

| Module | Description |
|--------|-------------|
| `POPMIX_EXTRACTF2` | Extract F2 blocks from genotype data |
| `POPMIX_QPADM` | Estimate admixture proportions |
| `POPMIX_QPWAVE` | Test ancestry stream count |
| `POPMIX_DATES` | Estimate admixture timing |
| `POPMIX_CONVERTF` | Convert between genotype formats |
| `POPMIX_F2` | Compute f2 statistic |
| `POPMIX_F3` | Compute f3 statistic |
| `POPMIX_F4` | Compute f4 statistic |

## Quick Start

### 1. Include modules in your pipeline

```nextflow
include { POPMIX_EXTRACTF2 } from './popmix/extractf2/main'
include { POPMIX_QPADM     } from './popmix/qpadm/main'
include { POPMIX_DATES     } from './popmix/dates/main'
```

### 2. Define input channel

```nextflow
// Input: tuple of meta map and genotype files
ch_input = Channel.of([
    [ id: 'sample1', prefix: '/path/to/data' ],
    file('/path/to/data.*')
])
```

### 3. Run modules

```nextflow
// Extract F2 blocks
POPMIX_EXTRACTF2(
    ch_input,
    ['French', 'Yoruba', 'Han'],  // populations
    0.05                           // block size
)

// Run qpAdm
POPMIX_QPADM(
    ch_input,
    'French',                      // target
    ['Yamnaya', 'Anatolia_N'],     // sources
    ['Mbuti', 'Han', 'Papuan']     // outgroups
)
```

## Input Format

The input channel should contain tuples of:
1. `meta` - A map with at least:
   - `id`: Sample identifier
   - `prefix`: File path prefix (without extension)
   - `format`: (optional) Input format hint
2. `geno_files` - Path to genotype files

### Example Samplesheet

```csv
sample,prefix,format
test1,/data/genotypes1,packedancestrymap
test2,/data/genotypes2,plink
```

### Supported Formats

- `plink` / `packedped`: Binary PLINK (.bed, .bim, .fam)
- `eigenstrat`: Text EIGENSTRAT (.geno, .snp, .ind)
- `packedancestrymap`: Binary EIGENSTRAT (.geno, .snp, .ind)
- `ped`: Text PLINK (.ped, .map)
- `ancestrymap`: Sparse format (.ancestrymapgeno, .snp, .ind)

## Output

Each module produces:
- Results in JSON format (for programmatic parsing)
- `versions.yml` for pipeline provenance

## Configuration

See `nextflow.config.example` for example configuration.

### Resource Labels

| Label | Description | Default |
|-------|-------------|---------|
| `process_low` | Light tasks (f2, f3, f4) | 2 CPU, 12 GB |
| `process_medium` | Medium tasks (extract-f2, dates) | 4 CPU, 36 GB |
| `process_high` | Heavy tasks (qpadm rotation) | 8 CPU, 72 GB |

## Example Workflow

See `example_workflow.nf` for a complete example pipeline.

```bash
nextflow run example_workflow.nf \
    --input samplesheet.csv \
    --populations "French,Yoruba,Han,Mbuti" \
    --target French \
    --sources "Yamnaya,Anatolia_N" \
    --outgroups "Mbuti,Han,Papuan"
```

## Container Images

Modules support both Docker and Singularity:

```nextflow
profiles {
    docker {
        docker.enabled = true
        process.container = 'biocontainers/popmix:0.1.0'
    }
    singularity {
        singularity.enabled = true
        process.container = 'docker://biocontainers/popmix:0.1.0'
    }
}
```

## Integration with nf-core Pipelines

These modules follow nf-core conventions and can be integrated into nf-core pipelines:

1. Copy modules to `modules/local/popmix/`
2. Update `modules.config` with process-specific options
3. Include modules in your workflow

## License

MIT License
