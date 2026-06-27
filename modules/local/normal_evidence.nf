// Pipeline A glue: pile up the matched buffy-coat BAM at the somatic sites to
// get per-site (normal_alt, normal_depth). panel_select.py uses this to subtract
// CHIP / germline (design P4 — CHIP is defined by presence in the WBC compartment).
//
// TODO: bin/normal_evidence.py is NOT written yet. The stub keeps the DAG valid.
// Real impl: bcftools mpileup -R <sites> the buffy BAM, emit the per-site table.
process NORMAL_EVIDENCE {
    tag "${meta.id}"
    label 'process_low'
    container 'mrd-umi/utils:1.0'

    input:
    tuple val(meta), path(somatic_vcf), path(vcf_tbi), path(buffy_bam), path(buffy_bai)
    path fasta

    output:
    tuple val(meta), path("*.normal_evidence.tsv"), emit: evidence
    path "versions.yml", emit: versions

    script:
    """
    normal_evidence.py pileup \\
        --vcf ${somatic_vcf} --bam ${buffy_bam} --ref ${fasta} \\
        --out ${meta.id}.normal_evidence.tsv
    echo '"${task.process}": {normal_evidence: TODO}' > versions.yml
    """

    stub:
    "touch ${meta.id}.normal_evidence.tsv versions.yml"
}
