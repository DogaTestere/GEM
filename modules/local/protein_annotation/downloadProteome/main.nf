process DOWNLOAD_PROTEOME_NCBI {
    tag "${meta.id}"
    label 'process_lowest'

    conda "bioconda::ncbi-datasets-cli"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-datasets-cli:16.15.0--pyhdfd78af_0' :
        'biocontainers/ncbi-datasets-cli:16.15.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.faa"), emit: pep
    tuple val("${task.process}"), val('datasets'), eval("datasets --version"), topic: versions, emit: versions_datasets
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    datasets download genome accession ${accession} \\
        --include protein \\
        --filename ${meta.id}_dataset.zip

    unzip -o ${meta.id}_dataset.zip
    find ncbi_dataset -name "protein.faa" -exec mv {} ${meta.id}.faa \\;
    """
}