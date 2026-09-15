/*
========================================================================================
    IMPORT NF-CORE MODULES
========================================================================================
*/

//
// MODULE: Local
//
include { FLAIR_COMBINE    } from '../../../modules/local/flair/combine/main'
include { FLAIR_QUANTIFY   } from '../../../modules/local/flair/quantify/main'




workflow FLAIR_COMBINE_QUANTIFY {
    
    take:
    flair_bed_ch
    flair_fa_ch
    flair_map_ch
    flair_gtf_ch
    filtered_reads_ch // filtered long reads, fed into QUANTIFY


    main:
    
    /*
     * Key all per-sample FLAIR outputs by sample id.
     */
    flair_bed_keyed = flair_bed_ch
        .map { meta, bed -> tuple(meta.id, meta, bed) }

    flair_fa_keyed = flair_fa_ch
        .map { meta, fa -> tuple(meta.id, fa) }

    flair_map_keyed = flair_map_ch
        .map { meta, map -> tuple(meta.id, map) }

    flair_gtf_keyed = flair_gtf_ch
        .map { meta, gtf -> tuple(meta.id, gtf) }


    /*
     * Join per-sample FLAIR outputs safely by sample id.
     */
    flair_samples = flair_bed_keyed
        .join(flair_fa_keyed)
        .join(flair_map_keyed)
        .join(flair_gtf_keyed)
        .map { sample_id, meta, bed, fa, map, gtf ->
            tuple(meta.group_id, meta, bed, fa, map, gtf)
        }

    /*
     * Group samples by group_id.
     */
    flair_by_group = flair_samples
        .groupTuple(by: 0)
        .map { group_id, metas, beds, fas, maps, gtfs ->
            def group_meta = [
                id      : group_id,
                group_id: group_id
            ]

            tuple(group_meta, metas, beds, fas, maps, gtfs)
        }


    flair_multi  = flair_by_group.filter { meta, metas, beds, fas, maps, gtfs -> metas.size() >= 2 }
    flair_single = flair_by_group.filter { meta, metas, beds, fas, maps, gtfs -> metas.size() == 1 }



    flair_multi_for_combine = flair_multi
    .map { meta, metas, beds, fas, maps, gtfs ->
        tuple(meta, metas, beds, fas, maps)
    }

    /*
     * Multi-sample groups: combine FLAIR isoforms.
     */
    FLAIR_COMBINE(flair_multi_for_combine)



    /*
     * Singleton groups: use original per-sample FLAIR outputs as combined outputs.
     */
    single_combined_gtf = flair_single.map { meta, metas, beds, fas, maps, gtfs ->
        tuple(meta, gtfs[0])
    }

    single_combined_fa = flair_single.map { meta, metas, beds, fas, maps, gtfs ->
        tuple(meta, fas[0])
    }

    single_combined_bed = flair_single.map { meta, metas, beds, fas, maps, gtfs ->
        tuple(meta, beds[0])
    }


    /*
     * Combined transcriptomes from both paths.
     */
    combined_gtf_ch = FLAIR_COMBINE.out.gtf.mix(single_combined_gtf)
    combined_fa_ch  = FLAIR_COMBINE.out.fa.mix(single_combined_fa)
    combined_bed_ch = FLAIR_COMBINE.out.bed.mix(single_combined_bed)



    /*
     * Group filtered long reads by group_id for quantification.
     */
    long_reads_by_group = filtered_reads_ch
        .map { meta, fastq ->
            tuple(meta.group_id, meta, fastq)
        }
        .groupTuple(by: 0)


    /*
     * Quantify both multi-sample and singleton groups.
     */
    flair_quantify_input = combined_fa_ch
        .map { meta, comb_fa ->
            tuple(meta.group_id, meta, comb_fa)
        }
        .join(long_reads_by_group)
        .map { group_id, group_meta, comb_fa, read_metas, fastqs ->
            tuple(group_meta, read_metas, fastqs, comb_fa)
        }

    FLAIR_QUANTIFY(flair_quantify_input)



    emit:
    quantify_counts = FLAIR_QUANTIFY.out.counts

    combined_gtf = combined_gtf_ch
    combined_fa  = combined_fa_ch
    combined_bed = combined_bed_ch
}