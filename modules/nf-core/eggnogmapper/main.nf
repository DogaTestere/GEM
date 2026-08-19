process EGGNOGMAPPER {
    tag "$meta.id"
    label 'process_high'

    // Dynamic storeDir based on genome_type
    storeDir { 
        def genome_type = meta.genome_type ?: 'all'
        def dir = params.eggnog_data_dirs[genome_type]
        if (!dir) {
            log.warn "No eggnog data directory found for genome_type: ${genome_type}. Using 'all'."
            dir = params.eggnog_data_dirs.all
        }
        return dir
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2':
        'quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(search_mode), path(db)

    output:
    tuple val(meta), path("*.emapper.annotations")   , emit: annotations
    tuple val(meta), path("*.emapper.seed_orthologs"), emit: orthologs, optional: true
    tuple val(meta), path("*.emapper.hits")          , emit: hits     , optional: true
    tuple val("${task.process}"), val('eggnog-mapper'), eval("emapper.py --version 2>&1 | grep -o 'emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' | sed 's/emapper-//'"), topic: versions, emit: versions_eggnogmapper
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.extension == '.gz'
    def fasta_name    = is_compressed ? fasta.baseName : "${fasta}"
    
    // Map search modes to their flags
    def db_flags = [
        'diamond': '--dmnd_db',
        'novel_fams': '--dmnd_db',
        'mmseqs': '--mmseqs_db',
        'hmmer': '--database',
        'no_search': '--annotate_hits_table',
        'cache': '--cache'
    ]
    
    // Handle db path
    def db_path = db ?: ''
    def db_arg = ''
    if (db && db_flags[search_mode]) {
        def db_path_formatted = db.isDirectory() ? db.toString() : db
        db_arg = "${db_flags[search_mode]} ${db_path_formatted}"
    }
    
    // Use --dbmem if memory > 40GB
    def dbmem = task.memory.toMega() > 40000 ? '--dbmem' : ''
    
    // For no_search and cache modes, db might be optional
    def skip_db_check = search_mode in ['no_search', 'cache'] && !db

    // Determine download options based on genome_type
    def db_type = meta.genome_type ?: 'all'
    def download_options = ''
    if (db_type == 'pro') {
        download_options = '--bacteria'
    } else if (db_type == 'euk') {
        download_options = '--eukaryota'
    } else {
        download_options = ''  // Download full database
    }

    // Get the storeDir path (which is the database directory)
    def data_dir = task.storeDir.toString()

    // Auto-download database if missing
    def downloadDb = """
    # Use a lock file to prevent concurrent downloads
    LOCK_FILE="${data_dir}/.download.lock"

    wait_for_lock() {
        while [ -f "\$LOCK_FILE" ]; do
            echo "Waiting for another download to complete..."
            sleep 10
        done
    }

    wait_for_lock
    touch "\$LOCK_FILE"

    if [ ! -f ${data_dir}/eggnog.db ]; then
        echo "EggNOG database not found at ${data_dir}. Downloading..."
    
        if ! command -v download_eggnog_data.py &> /dev/null; then
            echo "ERROR: download_eggnog_data.py not found. Please install eggnog-mapper."
            rm -f "\$LOCK_FILE"
            exit 1
        fi
    
        # Fix URL issue
        DOWNLOAD_SCRIPT=\$(which download_eggnog_data.py)
        if grep -q "eggnogdb.embl.de" \$DOWNLOAD_SCRIPT; then
            echo "Patching download script with updated URLs..."
            sed -i 's|eggnogdb.embl.de|eggnog5.embl.de|g' \$DOWNLOAD_SCRIPT
        fi
    
        echo "Downloading database for genome_type: ${db_type}"
        download_eggnog_data.py \\
            --data_dir ${data_dir} \\
            ${download_options} \\
            -y \\
            --force
        
        if [ \$? -ne 0 ]; then
            echo "ERROR: Failed to download EggNOG database"
            rm -f "\$LOCK_FILE"
            exit 1
        fi
        echo "Database download completed."
    else
        echo "EggNOG database found at ${data_dir}. Skipping download."
    fi

    rm -f "\$LOCK_FILE"
    """

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    ${downloadDb}

    # Run eggnog-mapper
    emapper.py \\
        --cpu ${task.cpus} \\
        -i ${fasta_name} \\
        --data_dir ${data_dir} \\
        -m ${search_mode} \\
        ${db_arg} \\
        ${dbmem} \\
        ${args} \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version 2>&1 | grep -o 'emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' | sed 's/emapper-//')
    END_VERSIONS
    """
}