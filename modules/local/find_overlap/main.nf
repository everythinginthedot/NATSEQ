process FIND_OVERLAP {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"

    input:
    // gtf + sqanti_classification MUST travel together in one tuple, keyed by the
    // same meta — they come from two independent upstream paths of different
    // length (classification straight from SQANTI3_QC; gtf only after
    // FILTER -> RESCUE -> AGAT_CONVERTSPGXF2GXF) that don't necessarily emit
    // same-sample items in the same relative order. Two separate positional
    // inputs would pair them by arrival order, not by meta.id (see the same bug,
    // same fix, in modules/local/sqanti3/qc/main.nf).
    tuple val(meta), path(gtf), path(sqanti_classification)
    path(ref_gtf)

    output:
    tuple val(meta), path("${meta.id}.nat.overlaps.tsv"), emit: main_tsv
    tuple val(meta), path("${meta.id}.nat.REF.overlaps.tsv"), emit: ref_tsv
    tuple val(meta), path("${meta.id}.nat_caller.log"), emit: nat_log


    script:
    """
    nat_caller.py \
        --query ${gtf} \
        --reference ${ref_gtf} \
        --output ${meta.id}.nat.overlaps.tsv \
        --reference-hints-output ${meta.id}.nat.REF.overlaps.tsv \
        --sqanti-classification ${sqanti_classification}
        > ${meta.id}.nat_caller.log 2>&1
    """
}