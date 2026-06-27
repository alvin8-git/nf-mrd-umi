// B3: build duplex consensus reads (require both strands to agree). For a
// single-strand assay swap this module for CallMolecularConsensusReads.
process FGBIO_CALLDUPLEXCONSENSUS {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/fgbio:4.1.0--hdfd78af_0'

    input:
    tuple val(meta), path(grouped_bam)

    output:
    tuple val(meta), path("*.consensus.unmapped.bam"), emit: bam
    path "versions.yml",                               emit: versions

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g CallDuplexConsensusReads \\
        --input ${grouped_bam} \\
        --min-reads ${params.consensus_min_reads} \\
        --threads ${task.cpus} \\
        --output ${meta.id}.consensus.unmapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | sed 's/^.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.consensus.unmapped.bam versions.yml"
}
