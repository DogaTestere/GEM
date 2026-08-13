process TRANSPORT_ADDER {
        tag "${meta.id}"
    label 'process_lowest'

    conda "$ {moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/library/python:3.11.9' }"

    input:
    tuple val(meta), path(initial_model), path(trans_met_db)
    path (python_script)

    output:
    tuple val(meta), path("*_transport_added_model"), emit:trans_added_model
    path "versions.yml", emit:versions

    script:
    """
    python3 ${python_script} \
        --input_db ${trans_met_db} \
        --input_model ${initial_model} \
        --output_model ${meta.id}_transport_added_model.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS 
    """
}