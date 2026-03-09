include { UNIPROT_MAPPING } from "../../../modules/local/uniprot"
include { GO_TERM_FINDER  } from "../../../modules/local/gene_onthology"
include { KEGG_REQUESTS   } from "../../../modules/local/kegg"
include { JSON_MERGING    } from "../../../modules/local/json_merging"

workflow WEB_REQUESTS {
  take:
    gbk_ch // tuple val(meta), path(gbk_json)

  main:
    mapping_script = channel.fromPath(params.mapping_script)
    UNIPROT_MAPPING(gbk_ch, mapping_script) 

    goTerm_script = channel.fromPath(params.goTerm_script)
    GO_TERM_FINDER(gbk_ch, goTerm_script)

    merged_kegg = gbk_ch
        .join(UNIPROT_MAPPING.out.mapping)

    kegg_requests_script = channel.fromPath(params.kegg_requests_script)
    KEGG_REQUESTS(merged_kegg, kegg_requests_script)
    
    merged_input = gbk_ch
        .join(GO_TERM_FINDER.out.go_terms)
        .join(KEGG_REQUESTS.out.kegg_db)
        
    channel.fromPath(params.go_basic).set { go_basic_file }
        
    json_merging_script = channel.fromPath(params.json_merging_script)
    JSON_MERGING(merged_input, json_merging_script, go_basic_file)

  emit:
    merged_db = JSON_MERGING.out.merged_db
    versions    = JSON_MERGING.out.versions
}