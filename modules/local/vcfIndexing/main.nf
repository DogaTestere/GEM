process BCFTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_lowest'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/kegg-requests:3.11' }"
        
    input:
    tuple val(meta), path(variant_bcf)
    
    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: indexed_vcf
    
    script:
    """
    bcftools convert -O v -o ${meta.id}.vcf ${variant_bcf}
    bgzip -c ${meta.id}.vcf > ${meta.id}.vcf.gz
    tabix -p vcf ${meta.id}.vcf.gz
    """
}