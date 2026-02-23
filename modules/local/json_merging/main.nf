process JSON_MERGING {
    tag "${meta.id}"
    label 'process_lowest'
    // irem bunu lowest yap sen

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/library/python:3.11.9' }"

    input:
    tuple val(meta), path(parsed_json), path(go_json), path(kegg_db)
    path (python_script)

    output:
    tuple val(meta), path("*_merged.db"), emit:merged_db
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --parsed_json ${parsed_json} \
        --go_json ${go_json} \
        --kegg_db ${kegg_db} \
        --added_db ${meta.id}_merged.db
    """
}