process BCFTOOLS_CALL {
    tag "${meta.id}"
    label 'process_lowest'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/kegg-requests:3.11' }"

    input:
    tuple val(meta), path(sorted_bam), path(indexed_bam), path(ref_file)
    
    output:
    tuple val(meta), path("*.bcf"), emit: variant_bcf
    
    script:
    """
    bcftools mpileup -f ${ref_file} ${sorted_bam} | bcftools call -mv -Ob -o ${meta.id}.bcf
    """
}