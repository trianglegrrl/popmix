process POPMIX_QPADM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popmix:0.1.0--py39h2e5f10c_0' :
        'biocontainers/popmix:0.1.0--py39h2e5f10c_0' }"

    input:
    tuple val(meta), path(geno_files)
    val(target)
    val(sources)
    val(outgroups)

    output:
    tuple val(meta), path("${prefix}_qpadm.json"), emit: results
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def left_arg = sources.join(',')
    def right_arg = outgroups.join(',')
    """
    popmix --threads ${task.cpus} qpadm \\
        --input ${meta.prefix} \\
        --target ${target} \\
        --left ${left_arg} \\
        --right ${right_arg} \\
        --format json \\
        ${args} > ${prefix}_qpadm.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"target":"${target}","sources":[],"weights":[]}' > ${prefix}_qpadm.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """
}
