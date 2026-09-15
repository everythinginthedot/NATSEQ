/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

//
// MODULE: Local
//

include { FLAIR } from '../../../modules/local/flair/transcriptome/main'  



workflow FLAIR_TRANSCRIPTOME {

    take:
    bam_ch
    junc_bed
    star_sj_by_group_ch
    genome_fa_ch
    annotation_ch


    main:

    // ========================================================================================
    // COMBINE LONG-READ AND OPTIONAL SHORT-READ JUNCTIONS FOR FLAIR
    // ========================================================================================

    /*
     * Long-read evidence is always present:
     *   bam_ch   : tuple(meta, bam, bai)
     *   junc_bed : tuple(meta, bed)
     *
     * After join:
     *   tuple(meta, bam, bai, bed)
     */
    bam_keyed = bam_ch
        .map { meta, bam, bai ->
            tuple(meta.id, meta, bam, bai)
        }

    junc_bed_keyed = junc_bed
        .map { meta, bed ->
            tuple(meta.id, bed)
        }

    long_input = bam_keyed
        .join(junc_bed_keyed)
        .map { sample_id, meta, bam, bai, bed ->
            tuple(meta, bam, bai, bed)
        }

    /*
     * Join by group_id. sj_map[group_id] gives this group's STAR junctions, or []
     * if none exist (e.g. --skip_aligning_star, or a tissue with no Illumina data).
     *
     * Uses combine() against a pre-collected map rather than
     * long_input.join(star_sj_by_group_ch, remainder: true): remainder:true buffers
     * the whole left channel until star_sj_by_group_ch itself closes, serializing
     * every sample behind the full IntronProspector batch instead of its own sample.
     */
    star_sj_map_ch = star_sj_by_group_ch
        .toList()
        .map { list -> list.collectEntries { group_id, sj -> [(group_id): sj] } }

    flair_input = long_input
        .map { meta, bam, bai, bed ->
            tuple(meta.group_id, meta, bam, bai, bed)
        }
        .combine(star_sj_map_ch)
        .map { group_id, meta, bam, bai, bed, sj_map ->
            tuple(meta, bam, bai, bed, sj_map[group_id] ?: [])
        }


    FLAIR(
        flair_input,
        genome_fa_ch,
        annotation_ch
    )


    emit:
    flair_fa      = FLAIR.out.fa
    flair_log     = FLAIR.out.log
    flair_counts  = FLAIR.out.counts
    flair_map     = FLAIR.out.map
    flair_bed     = FLAIR.out.bed
    flair_gtf     = FLAIR.out.gtf

}