process POPMIX_DATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popmix:0.1.0--py39h2e5f10c_0' :
        'biocontainers/popmix:0.1.0--py39h2e5f10c_0' }"

    input:
    tuple val(meta), path(geno_files)
    val(test_pop)
    val(source1)
    val(source2)

    output:
    tuple val(meta), path("${prefix}_dates.json"), emit: results
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    popmix --threads ${task.cpus} dates \\
        --input ${meta.prefix} \\
        --test-pop ${test_pop} \\
        --source1 ${source1} \\
        --source2 ${source2} \\
        --format json \\
        ${args} > ${prefix}_dates.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"test_pop":"${test_pop}","generations":0,"years_bp":0}' > ${prefix}_dates.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popmix: \$(popmix --version | sed 's/popmix, version //')
    END_VERSIONS
    """
}
