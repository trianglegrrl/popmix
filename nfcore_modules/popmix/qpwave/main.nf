process POPMIX_QPWAVE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popmix:0.1.0--py39h2e5f10c_0' :
        'biocontainers/popmix:0.1.0--py39h2e5f10c_0' }"

    input:
    tuple val(meta), path(f2_blocks)
    val(left_pops)
    val(right_pops)

    output:
    tuple val(meta), path("${prefix}_qpwave.json"), emit: results
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def left_str = left_pops.join(',')
    def right_str = right_pops.join(',')
    """
    popmix qpwave \\
        --f2-dir ${f2_blocks} \\
        --left ${left_str} \\
        --right ${right_str} \\
        --format json \\
        ${args} > ${prefix}_qpwave.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"rank":0,"p_value":1.0}' > ${prefix}_qpwave.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """
}
