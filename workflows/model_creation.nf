/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SPADES                 } from "../modules/nf-core/spades"
include { PROKKA                 } from '../modules/nf-core/prokka/main' 
include { MINIPROT_INDEX         } from '../modules/nf-core/miniprot/index/main' 
include { MINIPROT_ALIGN         } from '../modules/nf-core/miniprot/align/main'
include { PRODIGAL               } from '../modules/nf-core/prodigal/main' 
include { EGGNOGMAPPER           } from '../modules/nf-core/eggnogmapper/main'
include { GFFREAD                } from '../modules/nf-core/gffread/main'

// locale modules
include { DOWNLOAD_PROTEOME_NCBI } from '../modules/local/downloadProteome/'
include { GENEMARK_ES            } from "../modules/local/genemark_es"

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_model_creation_pipeline'

include { WEB_REQUESTS           } from "../subworkflows/local/web_requests"
include { MODEL_BUILDING         } from "../subworkflows/local/model_creation"
include { REFERENCE_ASSEMBLY     } from "../subworkflows/local/ref_assembly"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MODEL_CREATION {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_reads = ch_samplesheet
        .map { meta, fastqs -> tuple(meta, fastqs.collect { file(it) }) }

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    // FastQC 
    FASTQC(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions      = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Assembly
    //

    ch_reads
        .branch { meta, reads ->
            de_novo:    !meta.has_ref
            reference:   meta.has_ref
        }
        .set { ch_branched }

    // De-novo: SPAdes 
    SPADES(ch_branched.de_novo)

    // Reference assembly 
    REFERENCE_ASSEMBLY(
        ch_branched.reference.map { meta, reads ->
            tuple(meta, reads, meta.ref_file)   
        }
    )

    ch_contigs = SPADES.out.contigs.mix(REFERENCE_ASSEMBLY.out.contigs)

    //
    // Annotation
    //

    // This would need to further seperated if we change from the miniprot
    ch_contigs
        .branch { meta, reads -> 
            ab_initio: meta.ann_type == 'ab_in'
            reference: meta.ann_type == 'ref_in'
        }
        .set { ch_annotation}

    ch_annotation.ab_initio
        .branch { meta, reads ->
            pro: meta.genome_type == 'pro'
            euk: meta.genome_type == 'euk'
        }
        .set { ch_ab_initio }

    // Prokaryote ab-initio
    PRODIGAL(
        ch_ab_initio.pro,
        'gbk'
    )

    // Eukaryote ab-initio
    GENEMARK_ES(
        ch_ab_initio.euk,
        file(params.genemark_key)
    )

    // Turning genemark .gtf into .fasta
    ch_genemark_gff = GFFREAD(
        GENEMARK_ES.out.gtf.map { meta, gtf -> tuple(meta, gtf)},
        ch_ab_initio.euk.map { meta, fasta -> fasta}
    )

    // Proteome file checking and downloading
    ch_annotation.reference
        .branch { meta, reads ->
            has_pep : meta.pep_file
            down_pep : meta.pep_accession && !meta.pep_file
        }
        .set { ch_ref_split }

    // This is done so that later channel mix doesn't devolve into 5 line merge
    ch_pep_existing = ch_ref_split.has_pep_file
        .map { meta, reads ->
            tuple(meta, file(meta.pep_file))
        }

    DOWNLOAD_PROTEOME_NCBI(
        ch_ref_split.down_pep.map { meta, reads ->
            tuple(meta, meta.pep_accession)
        }
    )

    ch_complete_pep = ch_pep_existing.mix(DOWNLOAD_PROTEOME_NCBI.out.pep)

    // Miniprot
    MINIPROT_INDEX(
        ch_complete_pep
    )

    MINIPROT_ALIGN(
        ch_annotation.reference,
        MINIPROT_INDEX.out.index.map { meta, index -> tuple(meta, index) }
    )

    // Turning miniprot output into eggnogmapper readable version
    ch_gff = MINIPROT_ALIGN.out.gff
        .join(ch_annotation.reference)
    
    ch_miniprot_gff = GFFREAD(
        ch_gff.map { meta, gff, contigs -> tuple(meta, gff)},
        ch_gff.map { meta, gff, contigs -> contigs}
    )

    // Protein Functional Annotation
    ch_complete_annot = ch_genemark_gff.out.gffread_fasta
        .mix(ch_miniprot_gff.out.gffread_fasta)
        .mix(PRODIGAL.out.amino_acid_fasta)
    
    // This is for tuple val(search_mode), path(db)
    ch_eggnog_db = ch_complete_annot
        .map { meta, fasta ->
            def mode = meta.search_mode ?: 'diamond'
            def db_path = mode == 'diamond' ? params.eggnog_db_diamond :
                          mode == 'novel_fams' ? params.eggnog_db_novel_fams :
                          mode == 'mmseqs'  ? params.eggnog_db_mmseqs :
                          mode == 'hmmer'   ? params.eggnog_db_hmmer :
                          mode == 'no_search' ? params.eggnog_no_search_file :
                          params.eggnog_db_default
            tuple(mode, file(db_path))
        }
        .unique { it[0] }

    EGGNOGMAPPER(
        ch_complete_annot,
        ch_eggnog_db,
    )

    // 
    // WORKFLOW : Web Requests 
    //

    // !TODO: Adapt this to the eggnogmapper outputs

    WEB_REQUESTS(
        GB_PARSER.out.uni_first
    )

    //
    // WORKFLOW : Metabolic Model Creation
    //

    MODEL_BUILDING(
        WEB_REQUESTS.out.fixed_db
            .join(WEB_REQUESTS.out.go_terms)
    )

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'model_creation_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
