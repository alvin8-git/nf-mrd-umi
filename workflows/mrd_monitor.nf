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
    fasta_dict = file(params.fasta_dict)
    panel_vcf  = file(params.panel_vcf)
    panel_bed  = file(params.panel_bed)
    background = file(params.background)
    pon        = params.pon ? file(params.pon) : []

    UMI_CONSENSUS(ch_reads, fasta, fasta_dict, panel_bed)
    MRD_INTERROGATION(UMI_CONSENSUS.out.consensus_bam, panel_vcf, background, pon)

    MRD_INTERROGATION.out.report.view { meta, json -> "MRD ${meta.id} -> ${json}" }
}
