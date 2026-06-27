// B4 QC: on-target coverage of the collapsed consensus BAM. Converts the panel
// BED to a Picard interval_list inline (the picard image has BedToIntervalList).
process PICARD_COLLECTHSMETRICS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fasta_dict
    path  panel_bed

    output:
    tuple val(meta), path("*.hsmetrics.txt"), emit: metrics
    path "versions.yml",                      emit: versions

    script:
    """
    picard BedToIntervalList -I ${panel_bed} -SD ${fasta_dict} -O panel.interval_list
    picard CollectHsMetrics \\
        -I ${bam} -R ${fasta} \\
        -BI panel.interval_list -TI panel.interval_list \\
        -O ${meta.id}.hsmetrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectHsMetrics --version 2>&1 | head -1 | sed 's/^Version://')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.hsmetrics.txt versions.yml"
}
