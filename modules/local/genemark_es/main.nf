// modules/local/genemark_es.nf
process GENEMARK_ES {
    tag "$meta.id"
    label 'process_high'

    container 'ghcr.io/dogatestere/genmark-es:4.72'
    containerOptions = "-e HOME=/tmp/genemark_home"

    input:
    tuple val(meta), path(genome_fasta)
    path gm_key

    output:
    tuple val(meta), path("${prefix}/genemark.gtf"), emit: gtf
    tuple val(meta), path("${prefix}/output/gmhmm.mod"), emit: model, optional: true
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    mkdir -p \$HOME
    cp ${gm_key} \$HOME/.gm_key

    mkdir -p ${prefix}
    cd ${prefix}
    gmes_petap.pl \\
        --sequence ../${genome_fasta} \\
        --ES \\
        --cores ${task.cpus} \\
        ${args}
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genemark-es: 4.72
    END_VERSIONS
    """
}

