process MICROBEANNOTER_DOWNLOADER {
    tag "microbeannotator_db"
    label 'process_lowest'

    conda "${moduleDir}/environment.yml"
    container "..."

    // Cache location — set via params
    storeDir params.microbe_annoter_db

    input:
    val search_method   // 'diamond' (recommended) or 'blast' or 'sword'
    val use_light       // true to use --light

    output:
    path("db"), emit: db
    path "versions.yml", emit: versions

    script:
    def light_flag = use_light ? '--light' : ''
    """
    mkdir -p db
    microbeannotator_db_builder \\
        -d db \\
        -m ${search_method} \\
        -t ${task.cpus} \\
        ${light_flag} \\
        --no_aspera

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        microbeannotator: \$(microbeannotator --version 2>&1 | head -1)
    END_VERSIONS
    """
}