process SQANTI3_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path(classification)
    tuple val(meta3), path(ref_fasta)
    tuple val(meta4), path(ref_gtf)
    

    output:
    tuple val(meta), path("${prefix}_FILTER_output/"),          emit: results     // whole output dir
    tuple val(meta), path("${prefix}_FILTER_output/*.txt"),     emit: txt,        optional: true
    tuple val(meta), path("${prefix}_FILTER_output/*.pdf"),     emit: pdf,        optional: true
    tuple val(meta), path("${prefix}_FILTER_output/*.tsv"),     emit: tsv,        optional: true
    tuple val(meta), path("${prefix}_FILTER_output/*filtered.gtf"),     emit: gtf,        optional: true
    tuple val(meta), path("${prefix}_FILTER_output/*classification.txt"),     emit: classification,        optional: true
    path "versions.yml",                                        emit: versions


    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export PATH="${projectDir}/bin/SQANTI3:\$PATH"

    sqanti3_filter.py rules \
        ${args} \
        --sqanti_class ${classification} \
        --filter_gtf ${gtf} \
        --dir ${prefix}_FILTER_output \
        --output ${prefix}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_filter.py --version 2>&1 | sed 's/SQANTI3 version //')
    END_VERSIONS
    """
}
