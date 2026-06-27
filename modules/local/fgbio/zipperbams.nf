// Merge the aligned reads with the unmapped (u)BAM so the per-read tags
// (RX/UMI, or consensus tags) are restored onto the alignments. fgbio
// ZipperBams replaces the Picard MergeBamAlignment step.
process FGBIO_ZIPPERBAMS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fgbio:4.1.0--hdfd78af_0'

    input:
    tuple val(meta), path(aligned_sam), path(unmapped_bam)
    path  fasta            // reference FASTA (.dict alongside it)

    output:
    tuple val(meta), path("*.tagged.bam"), emit: bam
    path "versions.yml",                   emit: versions

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g ZipperBams \\
        --unmapped ${unmapped_bam} \\
        --input ${aligned_sam} \\
        --ref ${fasta} \\
        --output ${meta.id}.tagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | sed 's/^.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.tagged.bam versions.yml"
}
