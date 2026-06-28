// Pipeline B, stages B5-B6: interrogate the panel sites on the consensus BAM,
// then make the panel-integrated MRD call.
include { INTERROGATE   } from '../../modules/local/interrogate'
include { MRD_INTEGRATE } from '../../modules/local/mrd_integrate'

workflow MRD_INTERROGATION {
    take:
    ch_bam        // [ meta, bam, bai ]
    panel_vcf
    background
    pon           // [] if no panel-of-normals
    panel_lock    // [] if no patient-lock token

    main:
    INTERROGATE(ch_bam, panel_vcf)
    MRD_INTEGRATE(INTERROGATE.out.counts, background, pon, panel_vcf, panel_lock)

    emit:
    report = MRD_INTEGRATE.out.report   // [meta, *.mrd.json]
    counts = INTERROGATE.out.counts
}
