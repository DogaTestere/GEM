



workflow DATABASE_BUILDING {
  take:
    annot_ch // meta, pyrodigal/bakta/miniprot (maybe genemark)
    func_annot_ch // meta, metacerberus/microbeannoter/eggnog

  main:
    INITIAL_DB_CREATION() // takes the input channels and turns it into a .db

    TDB_2_SCRAPER() // takes a .db(can be any steps's db) and proteins.faa, returns a .db
    TCDB_SCRAPER() // takes the TDB_2 created .db, returns a .db

    GO_TERM_SEPERATOR() // takes a .db(can be any steps's db) and .obo file, returns a .db
    TRANSPORT_GUESSER() // takes seperator's .db, guesses substrates, returns a .db

    // !TODO: Find a way to run these based on the 'db_type'
    KEGG_BUILDER() // takes a .db(can be any steps's db) and gets reactions, returns a .db
    MODELSEED_BUILDER() // takes a .db(can be any steps's db) and gets reactions, returns a .db
    METACYC_BUILDER()// takes a .db(can be any steps's db) and gets reactions, returns a .db

    REACTION_BALANCING() // takes the last step's .db, returns a .db

  emit:
    // Balancer's output db
}