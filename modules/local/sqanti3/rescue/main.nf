process SQANTI3_RESCUE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(corrected_fasta), path(filtered_gtf), path(filt_class) // from QC + FILTER
    tuple val(meta2), path(ref_gtf)                                                 // annotation
    tuple val(meta3), path(ref_fasta)
    tuple val(meta4), path(ref_classification)                                      // from SQANTI3_QC_REF

    output:
    tuple val(meta), path("${prefix}_RESCUE_output/"), emit: results
    tuple val(meta), path("${prefix}_RESCUE_output/*_rescued.gtf"),            emit: gtf            
    tuple val(meta), path("${prefix}_RESCUE_output/*_rescued.fasta"),          emit: fasta,          optional: true
    path "versions.yml",                               emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    export PATH="${projectDir}/bin/SQANTI3:\$PATH"

    sqanti3_rescue.py rules \
        ${args} \
        --filter_class ${filt_class}\
        --rescue_isoforms ${corrected_fasta} \
        --rescue_gtf ${filtered_gtf} \
        --refGTF ${ref_gtf} \
        --refFasta ${ref_fasta} \
        --refClassif ${ref_classification} \
        --dir ${prefix}_RESCUE_output \
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_rescue.py --version 2>&1 | sed 's/SQANTI3 version //')
    END_VERSIONS
    """
}