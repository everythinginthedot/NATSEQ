/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

//
// MODULE: Local
//

include { SQANTI3_QC       } from '../../../modules/local/sqanti3/qc'  
include { SQANTI3_QC_REF   } from '../../../modules/local/sqanti3/qc_ref'  
include { SQANTI3_FILTER   } from '../../../modules/local/sqanti3/filter'  
include { SQANTI3_RESCUE   } from '../../../modules/local/sqanti3/rescue'  


workflow SQANTI {

    take:
    flair_cq_combined_gtf_ch
    genome_fa_ch
    annotation_ch
    flair_quantify_counts_ch   // FLAIR_QUANTIFY.out.counts, keyed by the same group meta as combined_gtf

    main:

    // gtf and counts MUST be joined into one channel (one tuple, one meta) before
    // being handed to SQANTI3_QC — passing them as two separate positional inputs
    // pairs items by arrival order, not by meta.id, and they come from two
    // independent upstream paths (combined_gtf vs quantify_counts) that don't
    // necessarily emit same-sample items in the same relative order.
    sqanti_input_gtf = flair_cq_combined_gtf_ch
        .map { meta, gtf -> tuple(meta.id, meta, gtf) }
        .join(flair_quantify_counts_ch.map { meta, counts -> tuple(meta.id, counts) })
        .map { group_id, meta, gtf, counts -> tuple(meta, gtf, counts) }


    SQANTI3_QC(
        sqanti_input_gtf,
        genome_fa_ch,
        annotation_ch
    )




    SQANTI3_FILTER(
        SQANTI3_QC.out.gtf,
        SQANTI3_QC.out.classification,
        genome_fa_ch,
        annotation_ch
    )




    
    
    if (params.annotation_classif == "") {

        SQANTI3_QC_REF(
            annotation_ch,
            genome_fa_ch
        )

        // ref_classification — meta-less value channel
        ref_classification = SQANTI3_QC_REF.out.classification.first()
    }

    else {

        ref_classification = channel.fromPath(params.annotation_classif).map { annotation_classif ->
            tuple([id: "ref_classification"], annotation_classif)
        }
        .first()

    }

    



    // --- SQANTI3 RESCUE ---
    // Join corrected_fasta (from QC) + filtered_gtf (from FILTER) by meta
    qc_fasta_keyed = SQANTI3_QC.out.fasta
        .map { meta, fasta ->
            tuple(meta.id, meta, fasta)
        }

    filter_gtf_keyed = SQANTI3_FILTER.out.gtf
        .map { meta, gtf ->
            tuple(meta.id, gtf)
        }

    filter_classification_keyed = SQANTI3_FILTER.out.classification
        .map { meta, classification ->
            tuple(meta.id, classification)
        }

    rescue_input = qc_fasta_keyed
        .join(filter_gtf_keyed)
        .join(filter_classification_keyed)
        .map { group_id, meta, fasta, gtf, classification ->
            tuple(meta, fasta, gtf, classification)
        }



    SQANTI3_RESCUE(
        rescue_input,
        annotation_ch,
        genome_fa_ch,
        ref_classification
    )

    
    emit:
    sqanti_qc_results      = SQANTI3_QC.out.results
    sqanti_classification  = SQANTI3_QC.out.classification
    sqanti_gtf_corrected   = SQANTI3_QC.out.gtf

    sqanti_filter_results = SQANTI3_FILTER.out.results

    sqanti_rescue_results = SQANTI3_RESCUE.out.results
    sqanti_gtf_rescued    = SQANTI3_RESCUE.out.gtf
}