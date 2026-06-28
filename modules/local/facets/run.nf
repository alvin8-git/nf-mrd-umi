// Pipeline A, A3: copy number + tumor purity/ploidy in one step (cnv_facets).
// Alternatives with the same role: PURPLE, TITAN (swap the module + container).
process FACETS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/cnv_facets:0.16.1--py312r43h8537716_1'

    input:
    tuple val(meta), path(tbam), path(tbai), path(nbam), path(nbai)
    path snp_vcf            // common SNPs (e.g. dbSNP) for the pileup
    path snp_tbi            // tabix index for snp_vcf (cnv_facets requires it alongside)

    output:
    tuple val(meta), path("*.facets.vcf.gz"), emit: cnv
    tuple val(meta), path("*.purity.txt"),    emit: purity
    path "versions.yml",                      emit: versions

    script:
    """
    cnv_facets.R --snp-tumour ${tbam} --snp-normal ${nbam} \\
        --snp-vcf ${snp_vcf} --out ${meta.id}.facets

    # purity/ploidy live in the output VCF header
    zcat ${meta.id}.facets.vcf.gz | grep -E '^##(purity|ploidy)=' > ${meta.id}.purity.txt || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnv_facets: 0.16.1
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.facets.vcf.gz ${meta.id}.purity.txt versions.yml"
}
