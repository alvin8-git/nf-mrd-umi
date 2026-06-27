// Pipeline A entry: matched tumor + buffy-coat WES -> personalized panel.
// Samplesheet CSV: patient,tumor_fastq_1,tumor_fastq_2,normal_fastq_1,normal_fastq_2
include { ALIGN_WES               } from '../subworkflows/local/align_wes'
include { GATK4_MUTECT2           } from '../modules/local/gatk4/mutect2'
include { GATK4_FILTERMUTECTCALLS } from '../modules/local/gatk4/filtermutectcalls'
include { FACETS                  } from '../modules/local/facets/run'
include { VEP                     } from '../modules/local/vep/annotate'
include { PYCLONE_PREP            } from '../modules/local/pyclone_prep'
include { PYCLONEVI               } from '../modules/local/pyclonevi/run'
include { NORMAL_EVIDENCE         } from '../modules/local/normal_evidence'
include { PANEL_SELECT            } from '../modules/local/panel_select'

workflow PANEL_DESIGN {
    // One row per patient -> two read groups (tumor + buffy), distinct sample names.
    ch_reads = Channel.fromPath(params.input)
        | splitCsv(header: true)
        | flatMap { row -> [
            tuple([ id: "${row.patient}_T", patient: row.patient, type: 'tumor',  sm: "${row.patient}_T" ],
                  [ file(row.tumor_fastq_1),  file(row.tumor_fastq_2)  ]),
            tuple([ id: "${row.patient}_N", patient: row.patient, type: 'normal', sm: "${row.patient}_N" ],
                  [ file(row.normal_fastq_1), file(row.normal_fastq_2) ])
        ] }

    fasta     = file(params.fasta)
    fai       = file(params.fasta_fai)
    dict      = file(params.fasta_dict)
    snp_vcf   = file(params.snp_vcf)
    chip_bl   = file(params.chip_blocklist)
    vep_cache = file(params.vep_cache)
    pon       = params.mutect2_pon       ? file(params.mutect2_pon)       : []
    germline  = params.germline_resource ? file(params.germline_resource) : []
    intervals = params.intervals         ? file(params.intervals)         : []

    ALIGN_WES(ch_reads, fasta)
    bams = ALIGN_WES.out.bam

    // pair tumor + normal by patient
    tumor  = bams.filter { it[0].type == 'tumor'  }.map { m, b, i -> [ m.patient, m, b, i ] }
    normal = bams.filter { it[0].type == 'normal' }.map { m, b, i -> [ m.patient, m, b, i ] }
    paired = tumor.join(normal).map { pt, tm, tb, ti, nm, nb, ni ->
        tuple([ id: pt, patient: pt, normal_sm: nm.sm ], tb, ti, nb, ni) }

    // A2: somatic discovery
    GATK4_MUTECT2(paired, fasta, fai, dict, pon, germline, intervals)
    GATK4_FILTERMUTECTCALLS(GATK4_MUTECT2.out.vcf, fasta, fai, dict)
    somatic = GATK4_FILTERMUTECTCALLS.out.vcf            // [meta, vcf, tbi]

    // A3: CN + purity/ploidy
    FACETS(paired, snp_vcf)

    // A5: annotation
    VEP(somatic, fasta, vep_cache)

    // A4: clonality / CCF  (prep is the remaining custom gap - see module TODO)
    ch_prep = somatic.join(FACETS.out.cnv).join(FACETS.out.purity)
        .map { m, vcf, tbi, cnv, pur -> tuple(m, vcf, tbi, cnv, pur) }
    PYCLONE_PREP(ch_prep)
    PYCLONEVI(PYCLONE_PREP.out.input)

    // buffy-coat evidence at the somatic sites (CHIP/germline subtraction)
    normal_for_ev = bams.filter { it[0].type == 'normal' }.map { m, b, i -> [ m.patient, b, i ] }
    ch_ne = somatic.map { m, vcf, tbi -> [ m.patient, m, vcf, tbi ] }
        .join(normal_for_ev)
        .map { pt, m, vcf, tbi, nb, ni -> tuple(m, vcf, tbi, nb, ni) }
    NORMAL_EVIDENCE(ch_ne, fasta)

    // A6: personalized panel
    ch_sel = VEP.out.vcf.join(PYCLONEVI.out.ccf).join(NORMAL_EVIDENCE.out.evidence)
        .map { m, vcf, tbi, ccf, ne -> tuple(m, vcf, tbi, ccf, ne) }
    PANEL_SELECT(ch_sel, chip_bl)

    PANEL_SELECT.out.panel.view { m, bed, vcf -> "PANEL ${m.id} -> ${bed}, ${vcf}" }
}
