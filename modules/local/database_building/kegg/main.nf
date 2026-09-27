process KEGG_REQUESTS {
    tag "${meta.id}"
    label 'process_lowest'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/kegg-requests:3.11' }"

    input:
    tuple val(meta), path(parsed_json), path(mapping_json)
    path (python_script)

    output:
    tuple val(meta), path("*.db"), emit:kegg_db
    path "versions.yml", emit:versions

    script:
    """
    python3 ${python_script} \
        --db ${meta.id}.db \
        --parsed_json ${parsed_json} \
        --kegg_map_json ${mapping_json} \
        --manual_json ${params.manual_fixes}
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}