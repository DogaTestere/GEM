include { INITIAL_MODEL   } from "../../../modules/local/init_model"

workflow MODEL_BUILDING {
    take:
        fixed_db_ch // tuple val(meta) path(fixed_db) path(go_terms)
       
    main:    
        initial_model_script = channel.fromPath(params.initial_model_script)
        INITIAL_MODEL(MODEL_BALANCING.out.balanced_db, initial_model_script)

        // Gapfilling belki

    emit:
        
}