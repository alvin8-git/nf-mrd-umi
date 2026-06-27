// B5: count unique consensus molecules supporting each panel mutation.
// Runs the custom engine (bin/interrogate.py) from the utils image.
process INTERROGATE {
    tag "${meta.id}"
    label 'process_medium'
    container 'mrd-umi/utils:1.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path  panel_vcf

    output:
    tuple val(meta), path("*.site_counts.tsv"), emit: counts
    path "versions.yml",                        emit: versions

    script:
    """
    interrogate.py run \\
        --bam ${bam} --panel ${panel_vcf} \\
        --min-bq ${params.min_bq} --min-mapq ${params.min_mapq} \\
        --out ${meta.id}.site_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interrogate: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.site_counts.tsv versions.yml"
}
