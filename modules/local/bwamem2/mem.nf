// Alignment. Emits an unsorted SAM (the bwa-mem2 biocontainer has no samtools,
// per the no-monolith rule); SAMTOOLS_SORT downstream makes the sorted BAM.
// -C copies the FASTQ comment (RX/UMI) onto the records; -p for interleaved.
process BWAMEM2_MEM {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_8'

    input:
    tuple val(meta), path(reads)
    path  fasta            // reference FASTA (bwa-mem2 index files alongside it)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml",            emit: versions

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA"
    """
    bwa-mem2 mem -t ${task.cpus} -p -C -R "${rg}" ${fasta} ${reads} > ${meta.id}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.sam versions.yml"
}
