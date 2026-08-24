process PSAMM_GENERATE_MODEL {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(eggnog_output)
    
    output:
        tuple val(meta), path("model.yaml"), emit: model
        tuple val(meta), path("compounds.yaml"), emit: compounds
        tuple val(meta), path("reactions.yaml"), emit: reactions

        // Nolur nolmaz, directory'de ouput olarak verelim
        tuple val(meta), path("${meta.id}_model"), emit: model_dir

    script:
    """
    psamm-generate-model generate-database \
        --annotations ${eggnog_output} \
        --type R \
        --out ${meta.id}_model
    """
}