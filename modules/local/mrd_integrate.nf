// B6: panel-integrated MRD call (empirical null + enrichment de-bias).
// Runs the custom engine (bin/mrd_integrate.py) from the utils image.
process MRD_INTEGRATE {
    tag "${meta.id}"
    label 'process_medium'
    container 'mrd-umi/utils:1.0'

    input:
    tuple val(meta), path(site_counts)
    path  background
    path  pon            // optional panel-of-normals (empirical null); [] if none
    path  panel_vcf

    output:
    tuple val(meta), path("*.mrd.json"), emit: report
    path "versions.yml",                 emit: versions

    script:
    def pon_arg = pon ? "--pon ${pon}" : ""
    """
    mrd_integrate.py run \\
        --site-counts ${site_counts} \\
        --background ${background} ${pon_arg} \\
        --panel ${panel_vcf} \\
        --patient-id ${meta.patient} --timepoint ${meta.timepoint} \\
        --min-molecules ${params.min_molecules} \\
        --out ${meta.id}.mrd.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mrd_integrate: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.mrd.json versions.yml"
}
