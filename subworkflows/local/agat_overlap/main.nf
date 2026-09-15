/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

//
// MODULE: Local
//

include { AGAT_CONVERTSPGXF2GXF } from '../../../modules/nf-core/agat/convertspgxf2gxf/main' 
include { FIND_OVERLAP          } from '../../../modules/local/find_overlap/main.nf' 


workflow AGAT_OVERLAP {
    take:
    sqanti_gtf_ch
    sqanti_classification_ch

    main:
    annotation = channel.fromPath(params.annotation).first()


    AGAT_CONVERTSPGXF2GXF (
        sqanti_gtf_ch
    )


    // Join by meta.id into one channel before FIND_OVERLAP — output_gff and
    // sqanti_classification_ch come from paths of different length (gtf went
    // through FILTER -> RESCUE -> AGAT; classification came straight from QC),
    // so two separate positional inputs would risk pairing by arrival order
    // instead of by sample.
    find_overlap_input = AGAT_CONVERTSPGXF2GXF.out.output_gff
        .map { meta, gtf -> tuple(meta.id, meta, gtf) }
        .join(sqanti_classification_ch.map { meta, classification -> tuple(meta.id, classification) })
        .map { group_id, meta, gtf, classification -> tuple(meta, gtf, classification) }


    FIND_OVERLAP (
        find_overlap_input,
        annotation
    )


    emit:

    rescued_gff = AGAT_CONVERTSPGXF2GXF.out.output_gff
    main_tsv    = FIND_OVERLAP.out.main_tsv
    ref_tsv     = FIND_OVERLAP.out.ref_tsv
    nat_log     = FIND_OVERLAP.out.nat_log
}