// modules/local/fastq_concat/main.nf
process FASTQ_CONCAT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(chunks)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: reads

    script:
    """
    cat ${chunks} > ${meta.id}_trimmed.fastq.gz
    """
}