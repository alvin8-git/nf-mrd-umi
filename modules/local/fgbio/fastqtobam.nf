// B1: FASTQ -> unmapped BAM with the UMI moved into the RX tag.
// read_structure is assay-specific: "8M+T 8M+T" for a duplex UMI design,
// "+T +T" for no UMI. Container assigned in conf/docker.config (FGBIO_.*).
process FGBIO_FASTQTOBAM {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fgbio:4.1.0--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.unmapped.bam"), emit: bam
    path "versions.yml",                     emit: versions

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g FastqToBam \\
        --input ${reads} \\
        --read-structures ${params.read_structure} \\
        --sample ${meta.id} --library ${meta.id} \\
        --output ${meta.id}.unmapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | sed 's/^.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.unmapped.bam versions.yml"
}
