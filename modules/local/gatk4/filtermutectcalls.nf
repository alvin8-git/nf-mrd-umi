// Pipeline A, A2: apply Mutect2's filters -> PASS somatic calls.
process GATK4_FILTERMUTECTCALLS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    gatk FilterMutectCalls -R ${fasta} -V ${vcf} --stats ${stats} \\
        -O ${meta.id}.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -i 'GATK' | head -1 | sed 's/.*v//')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.filtered.vcf.gz ${meta.id}.filtered.vcf.gz.tbi versions.yml"
}
