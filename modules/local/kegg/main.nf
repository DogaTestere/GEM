process KEGG_REQUESTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/kegg-requests:3.11' }"

    input:
    tuple val(meta), path(mapping_json)
    path (python_script)

    output:
    tuple val(meta), path("*.db"), emit:kegg_db
    path "versions.yml", emit:versions

    script:
    """
    python3 ${python_script} \
        --mapping_json ${mapping_json} \
        --db ${meta.id}.db
    """
}