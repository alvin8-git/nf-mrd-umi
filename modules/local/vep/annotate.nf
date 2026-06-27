// Pipeline A, A5: annotate somatic SNVs with Ensembl VEP (versioned cache ->
// auditable, regulatory-friendly). gnomAD population AF is added here too so
// panel_select can do population germline exclusion.
process VEP {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/ensembl-vep:116.0--pl5321h2a3209d_0'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path vep_cache

    output:
    tuple val(meta), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    vep --offline --cache --dir_cache ${vep_cache} \\
        --species homo_sapiens --assembly ${params.vep_assembly} \\
        --fasta ${fasta} --format vcf --vcf \\
        --compress_output bgzip \\
        -i ${vcf} -o ${meta.id}.vep.vcf.gz
    tabix -p vcf ${meta.id}.vep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl-vep: \$(vep --help 2>&1 | grep -i 'ensembl-vep' | head -1 | sed 's/.*: //')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.vep.vcf.gz ${meta.id}.vep.vcf.gz.tbi versions.yml"
}
