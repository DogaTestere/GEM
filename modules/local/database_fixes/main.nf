process MANUAL_CHECKS {
    tag "${meta.id}"
    label 'process_low'
    
    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/go-term-finder:3.11' }"
        
    input:
    tuple val(meta), path(kegg_db)
    path(python_script)
    
    output:
    tuple val(meta), path("*_manual.db"), emit: manual_fix
    path("versions.yml"), emit: versions
    
    script:
    """
    python3 ${python_script} \
        --kegg_db ${kegg_db} \
        --manual_inst ${params.manual_fixes} \
        --fixed_db ${meta.id}_manual.db
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}