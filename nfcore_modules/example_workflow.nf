#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Popmix Example Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Example Nextflow pipeline using popmix modules for population genetics analysis.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { POPMIX_EXTRACTF2 } from './popmix/extractf2/main'
include { POPMIX_F2        } from './popmix/f2/main'
include { POPMIX_F3        } from './popmix/f3/main'
include { POPMIX_F4        } from './popmix/f4/main'
include { POPMIX_QPADM     } from './popmix/qpadm/main'
include { POPMIX_QPWAVE    } from './popmix/qpwave/main'
include { POPMIX_DATES     } from './popmix/dates/main'
include { POPMIX_CONVERTF  } from './popmix/convertf/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // Example: Process genotype data
    //
    // Input channel: tuple of meta map and genotype files
    // meta must contain: id (sample ID), prefix (file prefix without extension)
    //
    // Example input samplesheet (CSV):
    // sample,prefix,format
    // test1,/path/to/data1,packedancestrymap
    // test2,/path/to/data2,plink

    ch_input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample,
                prefix: row.prefix,
                format: row.format ?: null
            ]
            def geno_files = file("${row.prefix}.*")
            [ meta, geno_files ]
        }

    // Define populations to analyze
    ch_populations = params.populations
        ? Channel.of(params.populations.split(',').toList())
        : Channel.of([])

    // Step 1: Extract F2 blocks
    POPMIX_EXTRACTF2(
        ch_input,
        ch_populations,
        params.block_size
    )

    // Step 2: Run qpAdm (if parameters provided)
    if (params.target && params.sources && params.outgroups) {
        ch_sources = params.sources.split(',').toList()
        ch_outgroups = params.outgroups.split(',').toList()

        POPMIX_QPADM(
            ch_input,
            params.target,
            ch_sources,
            ch_outgroups
        )
    }

    // Step 3: Run DATES (if parameters provided)
    if (params.test_pop && params.source1 && params.source2) {
        POPMIX_DATES(
            ch_input,
            params.test_pop,
            params.source1,
            params.source2
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Duration: $workflow.duration"
}
