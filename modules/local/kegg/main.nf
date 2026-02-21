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
    tuple val(meta), path("kegg.ready"), emit:kegg_finish
    path "versions.yml", emit:versions

    // Only file being send is a dummy file that makes later steps wait for this since sql databse is outside of work
    script:
    """
    mkdir -p "${params.kegg_cache}"

    python3 ${python_script} \
        --mapping_json ${mapping_json} \
        --db ${params.kegg_cache}/${meta.id}.db

    KEGG_CODE_HASH=\$(sha256sum ${python_script} | cut -d' ' -f1)
    echo "kegg_code=\${KEGG_CODE_HASH}" > kegg.ready
    """
}