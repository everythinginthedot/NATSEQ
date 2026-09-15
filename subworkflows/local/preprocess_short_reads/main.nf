/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//


include { FASTQC as FASTQC_RAW     } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../../../modules/nf-core/fastqc/main'
include { FASTP                    } from '../../../modules/nf-core/fastp/main'
include { TRIMMOMATIC              } from '../../../modules/nf-core/trimmomatic/main'



workflow PREPROCESS_SHORT_READS {
    take:
    short_reads_ch

    main:

    FASTQC_RAW( short_reads_ch )




    if (!params.skip_clipping) {
        if (params.clip_tool == 'fastp') {

            FASTP(
                short_reads_ch.map { meta, reads -> [meta, reads, []] },
                false,
                params.fastp_save_trimmed_fail,
                false
            )

            ch_short_reads_prepped = FASTP.out.reads

        }

        else if (params.clip_tool == 'trimmomatic') {

            TRIMMOMATIC( short_reads_ch )

            ch_short_reads_prepped = TRIMMOMATIC.out.trimmed_reads

        }
    }

    FASTQC_TRIMMED( ch_short_reads_prepped )



    emit:
    short_reads   = ch_short_reads_prepped

}
