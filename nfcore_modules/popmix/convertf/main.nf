process POPMIX_CONVERTF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popmix:0.1.0--py39h2e5f10c_0' :
        'biocontainers/popmix:0.1.0--py39h2e5f10c_0' }"

    input:
    tuple val(meta), path(geno_files)
    val(output_format)

    output:
    tuple val(meta), path("${prefix}.*"), emit: converted
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${output_format}"
    def input_format_arg = meta.format ? "-if ${meta.format}" : ''
    """
    popmix convertf \\
        --input ${meta.prefix} \\
        --output ${prefix} \\
        ${input_format_arg} \\
        --output-format ${output_format} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${output_format}"
    """
    touch ${prefix}.placeholder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """
}
