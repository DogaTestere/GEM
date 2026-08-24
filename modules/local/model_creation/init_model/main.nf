process INITIAL_MODEL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/library/python:3.11.9' }"

    input:
    tuple val(meta), path(balanced_db)
    path (python_script)

    output:
    tuple val(meta), path("*_initial.xml"), emit:initial_model
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --model_objective ${meta.rxn_id} \
        --model_direction ${meta.max_min} \
        --balanced_db ${balanced_db} \
        --output_model ${meta.id}_initial.xml
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}