// Pipeline B, stages B1-B4: raw reads -> molecularly-collapsed, realigned,
// coordinate-sorted consensus BAM (+ on-target QC). The alignment pattern
// (extract -> bwa -> zipper) runs twice, so those modules are imported under
// aliases (DSL2 requires a unique name per invocation).
//
// NOTE: fgbio ZipperBams expects the unmapped BAM queryname-grouped and the
// aligned input in matching order; the operational sort-order details still
// need validation against the fgbio best-practice flow (refinement, not wiring).

include { FGBIO_FASTQTOBAM                         } from '../../modules/local/fgbio/fastqtobam'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_RAW     } from '../../modules/local/samtools/fastq'
include { BWAMEM2_MEM    as BWAMEM2_RAW            } from '../../modules/local/bwamem2/mem'
include { FGBIO_ZIPPERBAMS as ZIPPER_RAW          } from '../../modules/local/fgbio/zipperbams'
include { FGBIO_GROUPREADSBYUMI                    } from '../../modules/local/fgbio/groupreadsbyumi'
include { FGBIO_CALLDUPLEXCONSENSUS                } from '../../modules/local/fgbio/callduplexconsensus'
include { FGBIO_FILTERCONSENSUSREADS              } from '../../modules/local/fgbio/filterconsensusreads'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_CONS    } from '../../modules/local/samtools/fastq'
include { BWAMEM2_MEM    as BWAMEM2_CONS           } from '../../modules/local/bwamem2/mem'
include { FGBIO_ZIPPERBAMS as ZIPPER_CONS          } from '../../modules/local/fgbio/zipperbams'
include { SAMTOOLS_SORT                            } from '../../modules/local/samtools/sort'
include { PICARD_COLLECTHSMETRICS                  } from '../../modules/local/picard/collecthsmetrics'

workflow UMI_CONSENSUS {
    take:
    ch_reads      // [ meta, [fastq_1, fastq_2] ]
    fasta
    fasta_dict
    panel_bed

    main:
    // B1: FASTQ -> uBAM (UMI -> RX)
    FGBIO_FASTQTOBAM(ch_reads)

    // B2: raw align (disposable; used only to assign UMI groups by coordinate)
    SAMTOOLS_FASTQ_RAW(FGBIO_FASTQTOBAM.out.bam)
    BWAMEM2_RAW(SAMTOOLS_FASTQ_RAW.out.reads, fasta)
    ch_zip_raw = BWAMEM2_RAW.out.sam.join(FGBIO_FASTQTOBAM.out.bam)   // [meta, sam, ubam]
    ZIPPER_RAW(ch_zip_raw, fasta)

    // B3: group by UMI + duplex consensus
    FGBIO_GROUPREADSBYUMI(ZIPPER_RAW.out.bam)
    FGBIO_CALLDUPLEXCONSENSUS(FGBIO_GROUPREADSBYUMI.out.bam)

    // B4: filter -> realign consensus -> sort
    FGBIO_FILTERCONSENSUSREADS(FGBIO_CALLDUPLEXCONSENSUS.out.bam, fasta)
    SAMTOOLS_FASTQ_CONS(FGBIO_FILTERCONSENSUSREADS.out.bam)
    BWAMEM2_CONS(SAMTOOLS_FASTQ_CONS.out.reads, fasta)
    ch_zip_cons = BWAMEM2_CONS.out.sam.join(FGBIO_FILTERCONSENSUSREADS.out.bam)
    ZIPPER_CONS(ch_zip_cons, fasta)
    SAMTOOLS_SORT(ZIPPER_CONS.out.bam)

    // B4: on-target QC
    PICARD_COLLECTHSMETRICS(SAMTOOLS_SORT.out.bam, fasta, fasta_dict, panel_bed)

    emit:
    consensus_bam = SAMTOOLS_SORT.out.bam               // [meta, bam, bai]
    hs_metrics    = PICARD_COLLECTHSMETRICS.out.metrics
}
