/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

//
// MODULE: Local
//

include { INTRONPROSPECTOR } from '../../../modules/local/intronprospector/main'  



workflow INTRONPROSPECTOR_JUNCTIONS {
    take:
    bam_ch
    genome_fa_unzip_ch
    genome_fai_ch

    main:

    bam_only = bam_ch.map { meta, bam, bai -> tuple(meta, bam) }

    INTRONPROSPECTOR(
        bam_only,      
        genome_fa_unzip_ch,
        genome_fai_ch
    )


    emit:
    junc_bed = INTRONPROSPECTOR.out.junctions

}