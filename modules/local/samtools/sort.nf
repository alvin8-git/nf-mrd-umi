// SAM -> coordinate-sorted, indexed BAM. Paired with BWAMEM2_MEM (whose
// biocontainer has no samtools, per the no-monolith rule).
process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    path "versions.yml",                           emit: versions

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.bam ${sam}
    samtools index -@ ${task.cpus} ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.bam ${meta.id}.bam.bai versions.yml"
}
