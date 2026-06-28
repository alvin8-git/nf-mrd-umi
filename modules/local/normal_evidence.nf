// Pipeline A glue: pile up the matched buffy-coat BAM at the somatic sites to
// get per-site (normal_alt, normal_depth). panel_select.py uses this to subtract
// CHIP / germline (design P4 — CHIP is defined by presence in the WBC compartment).
// Pileups the buffy BAM at the somatic sites -> per-site (normal_alt, normal_depth).
process NORMAL_EVIDENCE {
    tag "${meta.id}"
    label 'process_low'
    container 'mrd-umi/utils:1.0'

    input:
    tuple val(meta), path(somatic_vcf), path(vcf_tbi), path(buffy_bam), path(buffy_bai)
    path fasta
    path fai               // fasta .fai (pysam FastaFile needs it for --ref)

    output:
    tuple val(meta), path("*.normal_evidence.tsv"), emit: evidence
    path "versions.yml", emit: versions

    script:
    """
    normal_evidence.py pileup \\
        --vcf ${somatic_vcf} --bam ${buffy_bam} --ref ${fasta} \\
        --out ${meta.id}.normal_evidence.tsv
    echo '"${task.process}": {normal_evidence: bin}' > versions.yml
    """

    stub:
    "touch ${meta.id}.normal_evidence.tsv versions.yml"
}
