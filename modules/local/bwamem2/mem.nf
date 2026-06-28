// Alignment. Emits an unsorted SAM (the bwa-mem2 biocontainer has no samtools,
// per the no-monolith rule); SAMTOOLS_SORT downstream makes the sorted BAM.
// -C copies the FASTQ comment (RX/UMI) onto the records; -p for interleaved.
process BWAMEM2_MEM {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_8'

    input:
    tuple val(meta), path(reads)
    path  fasta            // reference FASTA
    path  index            // bwa-mem2 index sidecars (*.0123,*.amb,*.ann,*.bwt.2bit.64,*.pac), staged alongside fasta

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml",            emit: versions

    script:
    // -p for a single interleaved FASTQ (UMI flow); plain R1 R2 for WES.
    // -C carries the FASTQ comment (RX/UMI) when present; harmless otherwise.
    def sm = meta.sm ?: meta.id
    def rg = "@RG\\tID:${meta.id}\\tSM:${sm}\\tLB:${meta.id}\\tPL:ILLUMINA"
    def reads_list = reads instanceof List ? reads : [reads]
    def pe = reads_list.size() == 2 ? '' : '-p'
    """
    bwa-mem2 mem -t ${task.cpus} ${pe} -C -R "${rg}" ${fasta} ${reads_list.join(' ')} > ${meta.id}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.sam versions.yml"
}
