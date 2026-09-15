process SAMTOOLS_FAIDX {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/samtools:1.23--h96c455f_0"

    input:
    tuple val(meta), path(genome_fa)

    output:
    tuple val(meta), path("*.fai"), emit: fai
    tuple val(meta), path("*.{fa,fasta}"), emit: fa

    script:
    """
    gunzip -c ${genome_fa} > ${genome_fa.baseName}
    samtools faidx ${genome_fa.baseName}
    """
}