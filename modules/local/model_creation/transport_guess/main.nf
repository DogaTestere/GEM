process GUESS_TRANSPORT {
    tag "${meta.id}"
    label 'process_lowest'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/library/python:3.11.9' }"

    input:
    tuple val(meta), path(fixed_db), path (go_json)
    path (python_script)

    output:
    tuple val(meta), path("*_guessed.db"), emit:guess_db
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --input_db ${fixed_db} \
        --go_obo ${params.go_basic} \
        --go_json ${go_json} \
        --output_db ${meta.id}_guessed.db
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS 
    """
}