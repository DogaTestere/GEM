include { MODEL_BALANCING } from "../../../modules/local/model_balance"
include { INITIAL_MODEL   } from "../../../modules/local/init_model"

workflow MODEL_BUILDING {
    take:
        merged_db_ch // tuple val(meta) path(merged_db)
       
    main:
        //balancing_script = channel.fromPath(params.balancing_script)
        //MODEL_BALANCING(merged_db_ch, balancing_script)
        
        initial_model_script = channel.fromPath(params.initial_model_script)
        //INITIAL_MODEL(MODEL_BALANCING.out.balanced_db, initial_model_script)
        INITIAL_MODEL(merged_db_ch, initial_model_script)
        
        // Checklist modules if needed
        
    emit:
        finished_model = INITIAL_MODEL.out.initial_model
        versions = INITIAL_MODEL.out.versions
}