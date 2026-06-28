// Pipeline A, A6: select the personalized panel (SNV-only, CHIP + buffy +
// gnomAD germline filters, probe-feasibility-aware, truncal/CCF ranking).
// Runs the custom engine (bin/panel_select.py) from the utils image. Output is
// the BED + VCF that Pipeline B (mrd_monitor) interrogates.
process PANEL_SELECT {
    tag "${meta.id}"
    label 'process_low'
    container 'mrd-umi/utils:1.0'
    publishDir "${params.outdir}/panel_design", mode: 'copy'

    input:
    tuple val(meta), path(annotated_vcf), path(vcf_tbi), path(ccf_tsv), path(normal_evidence)
    path chip_blocklist

    output:
    tuple val(meta), path("*.panel.bed"), path("*.panel.vcf"), path("*.panel.lock"), emit: panel
    path "versions.yml", emit: versions

    script:
    """
    # PyClone-vi results -> per-variant ccf/clonal_prob panel_select expects
    pyclone_prep.py to-ccf --pyclone ${ccf_tsv} --out ccf.tsv

    panel_select.py run \\
        --vcf ${annotated_vcf} \\
        --ccf ccf.tsv \\
        --normal-evidence ${normal_evidence} \\
        --chip-blocklist ${chip_blocklist} \\
        --panel-size ${params.panel_size} --min-ccf ${params.min_ccf} \\
        --out-bed ${meta.id}.panel.bed --out-vcf ${meta.id}.panel.vcf

    # patient-lock: a content+patient hash of the panel; Pipeline B refuses to
    # interrogate a cfDNA sample against a panel whose lock does not match.
    sample_id.py provenance --panel ${meta.id}.panel.vcf \\
        --patient-id ${meta.id} --out ${meta.id}.panel.lock

    echo '"${task.process}": {panel_select: bin, sample_id: bin}' > versions.yml
    """

    stub:
    "touch ${meta.id}.panel.bed ${meta.id}.panel.vcf ${meta.id}.panel.lock versions.yml"
}
