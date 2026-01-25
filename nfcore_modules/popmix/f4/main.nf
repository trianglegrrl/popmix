process POPMIX_F4 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popmix:0.1.0--py39h2e5f10c_0' :
        'biocontainers/popmix:0.1.0--py39h2e5f10c_0' }"

    input:
    tuple val(meta), path(f2_blocks)
    val(pop_a)
    val(pop_b)
    val(pop_c)
    val(pop_d)

    output:
    tuple val(meta), path("${prefix}_f4.json"), emit: results
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    popmix f4 \\
        --f2-dir ${f2_blocks} \\
        --pop-a ${pop_a} \\
        --pop-b ${pop_b} \\
        --pop-c ${pop_c} \\
        --pop-d ${pop_d} \\
        --format json \\
        ${args} > ${prefix}_f4.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"pop_a":"${pop_a}","pop_b":"${pop_b}","pop_c":"${pop_c}","pop_d":"${pop_d}","f4":0,"se":0}' > ${prefix}_f4.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """
}
