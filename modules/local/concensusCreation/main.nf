process BCFTOOLS_CONCENSUS {
    tag "${meta.id}"
    label 'process_lowest'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/kegg-requests:3.11' }"
    
    input:
    tuple val(meta), path(vcf_file), path(vcf_index), path(ref_file)
    
    output:
    tuple val(meta), path("*_consensus.fasta"), emit:cons_fasta
    
    script:
    """
    bcftools consensus -f ${ref_file} ${vcf_file} > ${meta.id}_consensus.fasta
    """
}