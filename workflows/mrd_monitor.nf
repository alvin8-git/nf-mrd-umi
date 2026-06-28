// Pipeline B entry: cfDNA timepoint(s) -> MRD call.
// Samplesheet CSV columns: patient,timepoint,fastq_1,fastq_2
include { UMI_CONSENSUS     } from '../subworkflows/local/umi_consensus'
include { MRD_INTERROGATION } from '../subworkflows/local/mrd_interrogation'

workflow MRD_MONITOR {
    ch_reads = Channel.fromPath(params.input)
        | splitCsv(header: true)
        | map { row -> tuple(
            [ id: "${row.patient}_${row.timepoint}", patient: row.patient, timepoint: row.timepoint ],
            [ file(row.fastq_1), file(row.fastq_2) ]
          ) }

    fasta      = file(params.fasta)
    bwa_index  = files("${params.fasta}.{0123,amb,ann,bwt.2bit.64,pac}")
    fasta_dict = file(params.fasta_dict)
    panel_vcf  = file(params.panel_vcf)
    panel_bed  = file(params.panel_bed)
    background = file(params.background)
    pon        = params.pon ? file(params.pon) : []
    panel_lock = params.panel_lock ? file(params.panel_lock) : []

    UMI_CONSENSUS(ch_reads, fasta, fasta_dict, panel_bed, bwa_index)
    MRD_INTERROGATION(UMI_CONSENSUS.out.consensus_bam, panel_vcf, background, pon, panel_lock)

    MRD_INTERROGATION.out.report.view { meta, json -> "MRD ${meta.id} -> ${json}" }
}
