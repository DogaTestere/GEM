process METACERBERUS_ANNOTE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/metacerberus:1.4.0--pyhdfd78af_0'
        : 'biocontainers/metacerberus:1.4.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(faa_file)
    path(db) // Folder path from params. + /db

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    metacerberus.py \\
        --faa ${faa_file} \\
        --db_dir ${db} \\
        --outdir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacerberus: \$(metacerberus.py --version | sed 's/metacerberus //')
    END_VERSIONS
    """
}