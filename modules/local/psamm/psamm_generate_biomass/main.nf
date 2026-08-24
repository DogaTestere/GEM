process PSAMM_GENERATE_BIOMASS {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(model_dir)
        tuple val(meta2), path(genome_file)
        tuple val(meta3), path(proteome_file)
        tuple val(meta4), path(gff_file)
    
    output:
        tuple val(meta), path("model.yaml"), emit: biomass_model
        tuple val(meta), path("biomass_compounds.yaml"), emit: biomass_comps
        tuple val(meta), path("biomass_reactions.yaml"), emit: biomass_reacts

    script:
    """
    psamm-generate-model generate-biomass \
        --genome ${genome_file} \
        --proteome ${proteome_file} \
        --gff ${gff_file} \
        --model ${model_dir}
    """
}