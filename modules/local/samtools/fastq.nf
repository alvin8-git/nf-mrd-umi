// Extract reads from a (u)BAM as interleaved FASTQ, carrying the RX (UMI) tag
// into the read comment so bwa-mem2 -C can copy it back onto the alignments.
process SAMTOOLS_FASTQ {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.interleaved.fq.gz"), emit: reads
    path "versions.yml",                          emit: versions

    script:
    """
    samtools collate -@ ${task.cpus} -u -O ${bam} \\
        | samtools fastq -@ ${task.cpus} -T RX -0 /dev/null -s /dev/null - \\
        | gzip > ${meta.id}.interleaved.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.interleaved.fq.gz versions.yml"
}
