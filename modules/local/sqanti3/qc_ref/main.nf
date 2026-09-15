process SQANTI3_QC_REF {
    tag "ref"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(ref_gtf)
    tuple val(meta2), path(ref_fasta)

    output:
    tuple val(meta), path("ref_QC_output/"),                                         emit: results
    tuple val(meta), path("ref_QC_output/ref_classification.txt"),                   emit: classification, optional: true
    tuple val(meta), path("ref_QC_output/ref_corrected.gtf"),                        emit: gtf,            optional: true
    tuple val(meta), path("ref_QC_output/ref_corrected.fasta"),                      emit: fasta,          optional: true
    path "versions.yml",                                            emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    export PATH="${projectDir}/bin/SQANTI3:\$PATH"

    cp ${ref_gtf} annotation_isoforms.gtf

    sqanti3_qc.py \\
        ${args} \\
        --isoforms annotation_isoforms.gtf \\
        --refGTF ${ref_gtf} \\
        --refFasta ${ref_fasta} \\
        --dir ref_QC_output \\
        --output ref

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_qc.py --version 2>&1 | sed 's/SQANTI3 version //')
    END_VERSIONS
    """
}