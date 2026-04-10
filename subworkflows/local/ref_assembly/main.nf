include { BOWTIE_ALIGNMENT   } from "../../../modules/local/refAlignment"
include { BOWTIE_INDEXING    } from "../../../modules/local/refIndexing"
include { SAMTOOLS_SORT      } from "../../../modules/local/bamConversion"
include { SAMTOOLS_INDEX     } from "../../../modules/local/bamIndexing"
include { BCFTOOLS_CALL      } from "../../../modules/local/variantCalling"
include { BCFTOOLS_INDEX     } from "../../../modules/local/vcfIndexing"
include { BCFTOOLS_CONCENSUS } from "../../../modules/local/concensusCreation" 

workflow REFERENCE_ASSEMBLY {
    take:
        ref_channel // tuple (meta, reads, ref_file)

    main:
        ch_ref_only = ref_channel.map { meta, reads, ref_file ->
            tuple(meta, ref_file)
        }
        ch_reads_only = ref_channel.map { meta, reads, ref_file ->
            tuple(meta, reads)
        }

        BOWTIE_INDEXING(ch_ref_only)
        BOWTIE_ALIGNMENT(
            BOWTIE_INDEXING.out.indexed_ref
                .join(ch_reads_only)            // -> tuple(meta, index_dir, reads)
        )

        SAMTOOLS_SORT(
            BOWTIE_ALIGNMENT.out.aligned_bam,   
            ch_ref_only                          
        )
        SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sorted_bam)

        BCFTOOLS_CALL(
            SAMTOOLS_SORT.out.sorted_bam
                .join(SAMTOOLS_INDEX.out.indexed_bam)   // -> tuple(meta, bam, bai)
                .join(ch_ref_only)                      // -> tuple(meta, bam, bai, ref_file)
        )

        BCFTOOLS_INDEX(BCFTOOLS_CALL.out.variant_bcf)
        BCFTOOLS_CONCENSUS(
            BCFTOOLS_INDEX.out.indexed_vcf
                .join(ch_ref_only)                      // -> tuple(meta, vcf, vcf_index, ref_file)
        )

    emit:
        contigs = BCFTOOLS_CONCENSUS.out.cons_fasta  
}