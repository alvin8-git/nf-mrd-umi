// Pipeline A, A2: somatic discovery, tumor vs matched buffy-coat normal.
// PoN + germline resource are optional but recommended (filter artifacts/germline).
process GATK4_MUTECT2 {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'

    input:
    tuple val(meta), path(tbam), path(tbai), path(nbam), path(nbai)
    path fasta
    path fai
    path dict
    path pon            // [] if none
    path pon_tbi        // pon .tbi (GATK needs it alongside); [] if no pon
    path germline       // [] if none (e.g. gnomAD)
    path germline_tbi   // germline .tbi (GATK needs it alongside); [] if no germline
    path intervals      // [] for whole-exome default

    output:
    tuple val(meta), path("*.unfiltered.vcf.gz"), path("*.unfiltered.vcf.gz.tbi"), path("*.unfiltered.vcf.gz.stats"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def pon_arg  = pon       ? "--panel-of-normals ${pon}"   : ""
    def germ_arg = germline  ? "--germline-resource ${germline}" : ""
    def iv_arg   = intervals ? "-L ${intervals}"             : ""
    """
    gatk Mutect2 -R ${fasta} \\
        -I ${tbam} -I ${nbam} -normal ${meta.normal_sm} \\
        ${pon_arg} ${germ_arg} ${iv_arg} \\
        -O ${meta.id}.unfiltered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -i 'GATK' | head -1 | sed 's/.*v//')
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.unfiltered.vcf.gz ${meta.id}.unfiltered.vcf.gz.tbi ${meta.id}.unfiltered.vcf.gz.stats versions.yml"
}
