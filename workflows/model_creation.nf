/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_model_creation_pipeline'

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SPADES                 } from "../modules/nf-core/spades"
include { BAKTA_BAKTA            } from '../modules/nf-core/bakta/bakta/main'
include { BAKTA_BAKTADBDOWNLOAD  } from '../modules/nf-core/bakta/baktadbdownload/main'
include { MINIPROT_INDEX         } from '../modules/nf-core/miniprot/index/main' 
include { MINIPROT_ALIGN         } from '../modules/nf-core/miniprot/align/main'
include { PRODIGAL               } from '../modules/nf-core/prodigal/main' 
include { GFFREAD                } from '../modules/nf-core/gffread/main'

include { GUNZIP as GUNZIP_PRODIGAL_FAA } from '../modules/nf-core/gunzip/main' 
include { GUNZIP as GUNZIP_PRODIGAL_GFF } from '../modules/nf-core/gunzip/main' 

// locale modules
// protein annotation
include { METACERBERUS_DOWNLOAD     } from '../modules/local/protein_annotation/metacerberus/download'
include { METACERBERUS_ANNOTE       } from '../modules/local/protein_annotation/metacerberus/annote'
include { DOWNLOAD_PROTEOME_NCBI    } from '../modules/local/protein_annotation/downloadProteome'
include { MICROBEANNOTER_DOWNLOADER } from '../modules/local/protein_annotation/microbeannoter/download'
include { MICROBEANNOTER_ANNOTER    } from '../modules/local/protein_annotation/microbeannoter/annote'

// locale subworkflows
include { DATABASE_BUILDING      } from "../subworkflows/local/database_building"
include { MODEL_BUILDING         } from "../subworkflows/local/model_creation"
//include { PSAMM_MODEL_MAKING     } from "../subworkflows/local/psamm_model_making"

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

    // Kinda legacy? I want to keep the eggnogmapper as an option if we fix the rstudio issue
    ch_ab_initio.pro
        .branch { meta, reads ->
            eggnog: meta.func_ann == 'eggnog'
            default: meta.func_ann != 'eggnog'
        }
        .set { ch_ab_pro }

    //
    // Ab-initio default path

    if (params.bakta_db) {
        ch_bakta_db = Channel.value(file(params.bakta_db))
    } else {
        BAKTA_BAKTADBDOWNLOAD()
        ch_bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db.first()
    }

    BAKTA_BAKTA(
        ch_ab_pro.default,
        ch_bakta_db
    )

    METACERBERUS_DOWNLOAD(
        params.metacerberus_db_list
    )

    METACERBERUS_ANNOTE(
        BAKTA_BAKTA.out.faa,
        METACERBERUS_DOWNLOAD.out.db
    )

    // default ab-initio output channel
    ch_bakta_mc_ouput = BAKTA_BAKTA.out.faa
        .join(BAKTA_BAKTA.out.tsv)
        .join(METACERBERUS_ANNOTE.out.tsv)

    // 
    // Reference annotation default path

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

    MINIPROT_INDEX(
        ch_complete_pep
    )

    MINIPROT_ALIGN(
        ch_annotation.reference
            .join(MINIPROT_INDEX.out.index)
    )

    // Turns miniprot output to .fasta format so that it works both with microbeannoter and eggnogmapper
    ch_gff = MINIPROT_ALIGN.out.gff
        .join(ch_annotation.reference)
    
    ch_miniprot_gff = GFFREAD(
        ch_gff.map { meta, gff, contigs -> tuple(meta, gff)},
        ch_gff.map { meta, gff, contigs -> contigs}
    )

    MICROBEANNOTER_DOWNLOADER(
        params.microbe_annoter_search,
        params.microbe_annoter_light
    )

    ch_miniprot_gff.gffread_fasta
        .branch { meta, faa ->
            eggnog : meta.func_ann == 'eggnog'
            default: meta.func_ann != 'eggnog'
        }
        .set { ch_miniprot_by_annoter }
    
    MICROBEANNOTER_ANNOTER(
        ch_miniprot_by_annoter.default,
        MICROBEANNOTER_DOWNLOADER.out.db,
        params.microbe_annoter_search
    )

    // default ref_in ouput channel
    ch_mini_microbe_output = ch_miniprot_by_annoter.default
        .join(MICROBEANNOTER_ANNOTER.out.annotations)

    //
    // EGGNOG path

    PRODIGAL(
        ch_ab_pro.eggnog,
        'gff'
    )
    GUNZIP_PRODIGAL_FAA(PRODIGAL.out.amino_acid_fasta)
    // GUNZIP_PRODIGAL_GFF(PRODIGAL.out.gene_annotations)   // !TODO: re-enable if PSAMM needs it

    ch_eggnog_input = GUNZIP_PRODIGAL_FAA.out.gunzip
        .mix(ch_miniprot_by_annoter.eggnog)

    EGGNOGMAPPER(
        ch_eggnog_input,
        file(params.eggnog_db)
    )

    // eggnog output channel
    ch_eggnog_output = ch_eggnog_input
        .join(EGGNOGMAPPER.out.annotations)

    // 
    // WORKFLOW : Database building & ID conversions
    //

    DATABASE_BUILDING(
        ch_bakta_mc_ouput,
        ch_mini_microbe_output,
        ch_eggnog_output 
    )

    //
    // WORKFLOW : Metabolic Model Creation
    //

    MODEL_BUILDING(
        DATABASE_BUILDING.out.final_db
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
