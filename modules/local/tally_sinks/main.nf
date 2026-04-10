process TALLY_SINKS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/gb-parser:3.11' }"

    input:
    tuple val(meta), path(initial_model)
    path(python_script)

    output:
    tuple val(meta), path("*_tallied.xml"), emit: tallied_model
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --input_model ${initial_model} \
        --output_model "${meta.id}_tallied.xml" \
        --report_out "extra_report.txt"
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}