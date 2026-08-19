include { GUESS_TRANSPORT } from "../../../modules/local/transport_guess"
include { MODEL_BALANCING } from "../../../modules/local/model_balance"
include { INITIAL_MODEL   } from "../../../modules/local/init_model"

workflow MODEL_BUILDING {
    take:
        fixed_db_ch // tuple val(meta) path(fixed_db) path(go_terms)
       
    main:    
        // TODO: Adının değişmesi lazım, transport_selector gibi bişeyle
        tran_gues_script = channel.fromPath(params.tran_gues_script)
        GUESS_TRANSPORT(fixed_db_ch, tran_gues_script)
    
        // Tranport builder here

        // Transport adder here

        // Needs to take the output of transport adder
        balancing_script = channel.fromPath(params.balancing_script)
        MODEL_BALANCING(GUESS_TRANSPORT.out.guess_db, balancing_script)
        
        initial_model_script = channel.fromPath(params.initial_model_script)
        INITIAL_MODEL(MODEL_BALANCING.out.balanced_db, initial_model_script)

        // Gapfilling belki

    emit:
        
}