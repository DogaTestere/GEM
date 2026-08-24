process PSAMM_GENERATE_TRANSPORT {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(model_dir)
        tuple val(meta2), path(eggnog_output)
    
    output:
        tuple val(meta), path("model.yaml"), emit: transport_model
        tuple val(meta), path("${meta.id}_model"), emit: transport_dir
        tuple val(meta), path("transporters.yaml"), emit: transport_yaml

    script:
    """
    psamm-generate-model generate-transporters \
        --annotations ${eggnog_output} \
        --model ${model_dir}
    """
}