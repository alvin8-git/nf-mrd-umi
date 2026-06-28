// Pipeline A, A1: WES FASTQ pair -> aligned, duplicate-marked, sorted BAM.
// No UMIs here (unlike Pipeline B), so duplicates are marked, not collapsed.
include { BWAMEM2_MEM           } from '../../modules/local/bwamem2/mem'
include { SAMTOOLS_SORT         } from '../../modules/local/samtools/sort'
include { PICARD_MARKDUPLICATES } from '../../modules/local/picard/markduplicates'

workflow ALIGN_WES {
    take:
    ch_reads      // [ meta, [fastq_1, fastq_2] ]
    fasta
    bwa_index     // bwa-mem2 index sidecar files

    main:
    BWAMEM2_MEM(ch_reads, fasta, bwa_index)
    SAMTOOLS_SORT(BWAMEM2_MEM.out.sam)
    PICARD_MARKDUPLICATES(SAMTOOLS_SORT.out.bam)

    emit:
    bam = PICARD_MARKDUPLICATES.out.bam     // [meta, md.bam, md.bai]
}
