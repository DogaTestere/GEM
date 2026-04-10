process BOWTIE_ALIGNMENT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/go-term-finder:3.11' }"
        
    input:
        tuple val(meta), path(index_dir), path(fastq_reads)
    
    output:
        tuple val(meta), path("*.bam"), emit: aligned_bam // <-- Changed to .bam
    
    script:
    """
    bowtie2 \
        --threads ${task.cpus} \
        -x ${index_dir}/${meta.id} \
        -1 ${fastq_reads[0]} \
        -2 ${fastq_reads[1]} \
        | samtools view -bS - > ${meta.id}.bam
    """
}