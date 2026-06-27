// Pipeline A, A4: clonality / cancer-cell-fraction with PyClone-vi.
process PYCLONEVI {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/pyclone-vi:0.2.0--pyhdfd78af_0'

    input:
    tuple val(meta), path(pyclone_input)

    output:
    tuple val(meta), path("*.pyclone.tsv"), emit: ccf
    path "versions.yml", emit: versions

    script:
    """
    pyclone-vi fit -i ${pyclone_input} -o ${meta.id}.pyclone.h5 \\
        -c ${params.pyclone_clusters} -d beta-binomial
    pyclone-vi write-results-file -i ${meta.id}.pyclone.h5 -o ${meta.id}.pyclone.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyclone-vi: 0.2.0
    END_VERSIONS
    """

    stub:
    "touch ${meta.id}.pyclone.tsv versions.yml"
}
