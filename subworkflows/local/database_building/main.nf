include { BAKTA_INITIAL_DB_CREATION     } = from '../modules/local/database_building/bakta_db_creation'
include { MINIPROT_INITIAL_DB_CREATION  } = from '../modules/local/database_building/miniprot_db_creation'
include { EGGNOG_INITIAL_DB_CREATION    } = from '../modules/local/database_building/eggnog_db_creation'

include { KEGG_BUILDER                  } = from '../modules/local/database_building/kegg'
include { MODELSEED_BUILDER             } = from '../modules/local/database_building/modelseed'
include { METACYC_BUILDER               } = from '../modules/local/database_building/metacyc'

include { TDB_2_SCRAPER                 } = from '../modules/local/database_building/tdb2_scraper'
include { TCDB_SCRAPER                  } = from '../modules/local/database_building/tcdb_scraper'
include { GO_TERM_SEPERATOR             } = from '../modules/local/database_building/go_seperator'
include { TRANSPORT_GUESSER             } = from '../modules/local/database_building/transport_guesser'

include { REACTION_BALANCING            } = from '../modules/local/database_building/reaction_balancer'

workflow DATABASE_BUILDING {
  take:
    ch_bakta_mc_ouput // (meta, faa, bakta_tsv, mc_tsv)
    ch_mini_microbe_output // (meta, faa, ma_tsv)
    ch_eggnog_output // (meta, faa, eggnog_tsv)

  main:
    // Seperated bcs each onse has a different file structure
    // And imo it would be hard too much of a headache to make script accomadate everything instead of making three seperate ones
    BAKTA_INITIAL_DB_CREATION(ch_bakta_mc_ouput)
    MINIPROT_INITIAL_DB_CREATION(ch_mini_microbe_output) 
    EGGNOG_INITIAL_DB_CREATION(ch_eggnog_output)
    // Needs EC numbers, GO IDs, UniProt ID, Protein Name minimally
    // Other identifiers gotten from annotations are also needed

    // Any ID conversion step if needed

    // !TODO: Check these after writing the modules
    ch_initial_db = BAKTA_INITIAL_DB_CREATION.out.inital_db
      .mix(MINIPROT_INITIAL_DB_CREATION.out.initial_db)
      .mix(EGGNOG_INITIAL_DB_CREATION.out.initial_db)

    ch_initial_db
        .branch { meta, db, faa ->
            kegg:      meta.db_type == 'kegg'
            modelseed: meta.db_type == 'modelseed'
            metacyc:   meta.db_type == 'metacyc'
        }
        .set { ch_by_dbsource }
    
    // !TODO: Check if any ModelSEED path and MetaCyc would require something else
    KEGG_BUILDER(ch_by_dbsource.kegg)
    MODELSEED_BUILDER(ch_by_dbsource)
    METACYC_BUILDER(ch_by_dbsource)

    ch_reactions_db = KEGG_BUILDER.out.reactions_db
      .mix(MODELSEED_BUILDER.out.reactions_db)
      .mix(METACYC_BUILDER.out.reactions_db)

    TDB_2_SCRAPER(ch_reactions_db)  // (meta, db, faa) -> (meta, db) Bcs after this step we don't need .faa
    TCDB_SCRAPER(TDB_2_SCRAPER.out.tdb2_db)

    GO_TERM_SEPERATOR(TCDB_SCRAPER.out.tcdb_db, file(params.obo_file)) 
    TRANSPORT_GUESSER(GO_TERM_SEPERATOR.out.go_terms_db)

    REACTION_BALANCING(TRANSPORT_GUESSER.out.guessed_db) 

  emit:
    final_db = REACTION_BALANCING.out.final_db
}