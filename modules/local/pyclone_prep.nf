// Pipeline A glue: build the PyClone-vi input table from the somatic VCF +
// FACETS copy number + tumour purity, and map PyClone's output back to the
// per-variant (ccf, clonal_prob) TSV that panel_select.py expects. The `to-ccf`
// half runs inside PANEL_SELECT (both scripts live in the utils image).
process PYCLONE_PREP {
    tag "${meta.id}"
    label 'process_low'
    container 'mrd-umi/utils:1.0'

    input:
    tuple val(meta), path(somatic_vcf), path(vcf_tbi), path(facets_vcf), path(purity)

    output:
    tuple val(meta), path("*.pyclone_input.tsv"), emit: input
    path "versions.yml", emit: versions

    script:
    """
    pyclone_prep.py build \\
        --vcf ${somatic_vcf} --facets ${facets_vcf} --purity ${purity} \\
        --out ${meta.id}.pyclone_input.tsv
    echo '"${task.process}": {pyclone_prep: bin}' > versions.yml
    """

    stub:
    "touch ${meta.id}.pyclone_input.tsv versions.yml"
}
