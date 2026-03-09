process MODEL_BALANCING {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/library/python:3.11.9' }"

    input:
    tuple val(meta), path(merged_db)
    path (python_script)

    output:
    tuple val(meta), path("*_balanced.db"), emit:balanced_db
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --input ${merged_db} \
        --output ${meta.id}_balanced.db 
    """
}