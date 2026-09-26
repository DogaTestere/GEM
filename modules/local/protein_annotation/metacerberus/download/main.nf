process METACERBERUS_DOWNLOAD{
    tag "metacerberus_download"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/metacerberus:1.4.0--pyhdfd78af_0'
        : 'biocontainers/metacerberus:1.4.0--pyhdfd78af_0'}"

    storeDir params.metacerberus_db

    input:
        val databases

    output:
        path("db"), emit: db
        path "versions.yml", emit: versions

    script:
        def db_args = databases ? databases.join(' ') : ''
        """
        metacerberus.py --download ${db_args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            metacerberus: \$(metacerberus.py --version | sed 's/metacerberus //')
        END_VERSIONS
        """
}