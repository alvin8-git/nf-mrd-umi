// Pipeline A: mark PCR duplicates on WES (no UMIs here, unlike Pipeline B).
// --CREATE_INDEX emits the .bai so we don't need samtools in this image.
process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.md.bam"), path("*.md.bai"), emit: bam
    tuple val(meta), path("*.md.metrics"),               emit: metrics
    path "versions.yml",                                 emit: versions

    script:
    """
    picard -Xmx${task.memory.toGiga()}g MarkDuplicates \\
        -I ${bam} -O ${meta.id}.md.bam -M ${meta.id}.md.metrics --CREATE_INDEX true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | head -1 | sed 's/^Version://')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.md.bam ${meta.id}.md.bai ${meta.id}.md.metrics versions.yml"
}
