process GUESS_TRANS_COMPOUND {
    tag "${meta.id}"
    label 'process_lowest'

    conda "$ {moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/library/python:3.11.9' }"

    input:
    tuple val(meta), path(balanced_db)
    path (python_script)

    output:
    tuple val(meta), path("*_transport_met_guess.db"), emit:trans_met_db
    path "versions.yml", emit:versions

    script:
    """
    python3 ${python_script} \
        --input_db ${balanced_db}
        --output_db ${meta.id}_transport_met_guess.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS 
    """
}