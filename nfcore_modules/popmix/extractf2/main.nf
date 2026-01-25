process POPMIX_EXTRACTF2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popmix:0.1.0--py39h2e5f10c_0' :
        'biocontainers/popmix:0.1.0--py39h2e5f10c_0' }"

    input:
    tuple val(meta), path(geno_files)
    val(populations)
    val(block_size)

    output:
    tuple val(meta), path("${prefix}_f2"), emit: f2_blocks
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def pop_arg = populations ? "--populations ${populations.join(',')}" : ''
    def block_arg = block_size ? "--block-size ${block_size}" : ''
    """
    popmix --threads ${task.cpus} extract-f2 \\
        --input ${meta.prefix} \\
        --output ${prefix}_f2 \\
        ${pop_arg} \\
        ${block_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_f2
    touch ${prefix}_f2/metadata.json
    touch ${prefix}_f2/block_sums.npy
    touch ${prefix}_f2/block_weights.npy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """
}
