process POPMIX_F3 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popmix:0.1.0--py39h2e5f10c_0' :
        'biocontainers/popmix:0.1.0--py39h2e5f10c_0' }"

    input:
    tuple val(meta), path(f2_blocks)
    val(target)
    val(source1)
    val(source2)

    output:
    tuple val(meta), path("${prefix}_f3.json"), emit: results
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    popmix f3 \\
        --f2-dir ${f2_blocks} \\
        --target ${target} \\
        --source1 ${source1} \\
        --source2 ${source2} \\
        --format json \\
        ${args} > ${prefix}_f3.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"target":"${target}","source1":"${source1}","source2":"${source2}","f3":0,"se":0}' > ${prefix}_f3.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """
}
