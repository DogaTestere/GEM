process BOWTIE_INDEXING {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/go-term-finder:3.11' }"
    
    input:
    tuple val(meta), path(ref_file)

    output:
    tuple val(meta), path("ref_index"), emit: indexed_ref

    script:
    """
    mkdir -p ref_index
    bowtie2-build --threads ${task.cpus} ${ref_file} ref_index/${meta.id}
    """
}