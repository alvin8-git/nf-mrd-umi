// B4: strip consensus reads with weak support. min-reads is assay-specific
// (duplex takes three values: total per-strand-a per-strand-b).
process FGBIO_FILTERCONSENSUSREADS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fgbio:4.1.0--hdfd78af_0'

    input:
    tuple val(meta), path(consensus_bam)
    path  fasta

    output:
    tuple val(meta), path("*.filtered.bam"), emit: bam
    path "versions.yml",                     emit: versions

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g FilterConsensusReads \\
        --input ${consensus_bam} \\
        --ref ${fasta} \\
        --min-reads ${params.filter_min_reads} \\
        --min-base-quality ${params.filter_min_baseq} \\
        --output ${meta.id}.filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | sed 's/^.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.filtered.bam versions.yml"
}
