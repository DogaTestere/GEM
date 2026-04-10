include { GUESS_TRANSPORT } from "../../../modules/local/transport_guess"
include { MODEL_BALANCING } from "../../../modules/local/model_balance"
include { INITIAL_MODEL   } from "../../../modules/local/init_model"
include { TALLY_SINKS     } from "../../../modules/local/tally_sinks"

workflow MODEL_BUILDING {
    take:
        fixed_db_ch // tuple val(meta) path(fixed_db) path(go_terms)
       
    main:    
        tran_gues_script = channel.fromPath(params.tran_gues_script)
        GUESS_TRANSPORT(fixed_db_ch, tran_gues_script)
    
        balancing_script = channel.fromPath(params.balancing_script)
        MODEL_BALANCING(GUESS_TRANSPORT.out.guess_db, balancing_script)
        
        initial_model_script = channel.fromPath(params.initial_model_script)
        INITIAL_MODEL(MODEL_BALANCING.out.balanced_db, initial_model_script)
        
        //gapfilling_script = channel.fromPath(params.gapfilling_script)
        //GAPFILLING_STEP(INTIAL_MODEL.out.initial_model, gapfilling_script)
        
        tally_sinks_script = channel.fromPath(params.tally_sinks_script)
        TALLY_SINKS(INITIAL_MODEL.out.initial_model, tally_sinks_script)
        
    emit:
        finished_model = TALLY_SINKS.out.tallied_model
        versions = TALLY_SINKS.out.versions
}