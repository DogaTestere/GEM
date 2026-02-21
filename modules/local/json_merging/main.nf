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
    tuple val(meta), path("merge.ready"), emit:merge_finish
    path "versions.yml", emit: versions

    // Only file being emitted is a dummy file, model creation needs to wait for sql database update

    script:
    """
    python3 ${python_script} \
        --parsed_json ${parsed_json} \
        --go_json ${go_json} \
        --kegg_db ${params.kegg_cache}/${meta.id}.db

    MERGE_CODE_HASH=`sha256sum ${python_script} | cut -d" " -f1`
    echo "merge_code=\${MERGE_CODE_HASH}" > merge.ready
    """
}