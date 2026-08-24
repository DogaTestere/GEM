include { PSAMM_GENERATE_MODEL     } from "../../../modules/local/psamm/psamm_generate_model"
include { PSAMM_GENERATE_TRANSPORT } from "../../../modules/local/psamm/psamm_generate_transporters"
include { PSAMM_GENERATE_BIOMASS   } from "../../../modules/local/psamm/psamm_generate_biomass"

workflow PSAMM_MODEL_MAKING {
    take:
        eggnog_output_ch // tuple val(meta), path(eggnog.tsv)
        genome_ch // tuple val(meta), path(genome.fasta)
        proteome_ch // tuple val(meta), paath(proteome.fasta)
        gff_ch // tuple val(meta), path(annot.gff)

    // !QUESTION: Maybe merging them is better? At least for genome,proteome and gff since they will directly go to biomass
    // !TODO: İreme sormak lazım

    // PSAMM modifies the model on the run and wants a DIRECTORY as the output
    main:
        PSAMM_GENERATE_MODEL(eggnog_output_ch)

        PSAMM_GENERATE_TRANSPORT(PSAMM_GENERATE_MODEL.out.model_dir, eggnog_output_ch)

        PSAMM_GENERATE_BIOMASS(
            PSAMM_GENERATE_TRANSPORT.out.transport_dir,
            genome_ch,
            proteome_ch,
            gff_ch
        )
    
    emit:
        final_model = PSAMM_GENERATE_BIOMASS.out.biomass_model

        // nolur nolmaz
        biomass_compounds = PSAMM_GENERATE_BIOMASS.out.biomass_comps
        biomass_reactions = PSAMM_GENERATE_BIOMASS.out.biomass_reacts
        transport_yaml = PSAMM_GENERATE_TRANSPORT.out.transport_yaml
        model_compounds = PSAMM_GENERATE_MODEL.out.compounds
        model_reactions = PSAMM_GENERATE_MODEL.out.reactions
}