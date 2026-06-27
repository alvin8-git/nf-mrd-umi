// B3: group reads by UMI + mapped position. --strategy paired enables the
// duplex pairing (top/bottom strand); use 'adjacency' for single-strand (SSCS).
process FGBIO_GROUPREADSBYUMI {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fgbio:4.1.0--hdfd78af_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.grouped.bam"),    emit: bam
    tuple val(meta), path("*.grouped.hist.txt"), emit: histogram
    path "versions.yml",                       emit: versions

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g GroupReadsByUmi \\
        --input ${bam} \\
        --strategy ${params.umi_strategy} \\
        --output ${meta.id}.grouped.bam \\
        --family-size-histogram ${meta.id}.grouped.hist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | sed 's/^.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.grouped.bam ${meta.id}.grouped.hist.txt versions.yml"
}
