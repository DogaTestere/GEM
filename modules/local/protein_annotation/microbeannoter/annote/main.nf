process MICROBEANNOTER_ANNOTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/microbeannotator:2.0.5--pyhdfd78af_0'
        : 'biocontainers/microbeannotator:2.0.5--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(faa)
    path(db)
    val(search_method)

    output:
    tuple val(meta), path("annotation_results/*.annotations"), emit: annotations
    path "annotation_results/*.kofam"                       , emit: kofam_hits, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    microbeannotator \\
        -i ${faa} \\
        -d ${db} \\
        -o . \\
        -m ${search_method} \\
        -p 1 \\
        -t ${task.cpus} \\
        ${args}

    # Rename output to match nf-core conventions
    mkdir -p annotation_results
    mv ./*.annotations annotation_results/${prefix}.annotations 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        microbeannotator: \$(microbeannotator --version 2>&1 | sed 's/microbeannotator //')
    END_VERSIONS
    """
}