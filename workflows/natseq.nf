/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Installed directly from nf-core/modules
//
include { PREPROCESS_LONG_READS      } from '../subworkflows/local/preprocess_long_reads/main.nf'  
include { PREPROCESS_SHORT_READS     } from '../subworkflows/local/preprocess_short_reads/main.nf'                                                   
include { MINIMAP2                   } from '../subworkflows/local/minimap2/main.nf'       
include { FLAIR_TRANSCRIPTOME        } from '../subworkflows/local/flair_transcriptome/main.nf'    
include { INTRONPROSPECTOR_JUNCTIONS } from '../subworkflows/local/intronprospector_junctions/main.nf'   
include { FLAIR_COMBINE_QUANTIFY     } from '../subworkflows/local/flair_combine_quantify/main.nf'
include { SQANTI                     } from '../subworkflows/local/sqanti/main.nf'
include { AGAT_OVERLAP               } from '../subworkflows/local/agat_overlap/main.nf'
include { STAR_JUNCTIONS             } from '../subworkflows/local/intronprospector_star/main'


//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_FAIDX             } from '../modules/local/samtools/faidx/main'  
                                        



workflow NATSEQ {

    // Guard against a specific CLI-parsing footgun: `--flair_extra_args --stringent`
    // (a bare, argument-less flag with NO space, as its own shell token) gets
    // silently misparsed by Nextflow's own CLI parser as TWO separate params —
    // flair_extra_args becomes the literal string "true", and --stringent becomes
    // its own unused params.stringent — instead of flair_extra_args holding the
    // flag text. The run then silently succeeds using an unintended (and often
    // cache-hit-identical-to-another-broken-run) FLAIR config with no error at all.
    // Happened for real on 2026-08-25 (--stringent and --check_splice collided into
    // the same cache). Flags WITH a value (`--end_window 300`) are unaffected — the
    // embedded space stops Nextflow from mistaking them for a new flag. Correct
    // invocation either way: `--flair_extra_args=--stringent` (`=`, not a space).
    // (Lives here, not in nextflow.config: Nextflow's config parser v2 forbids
    // top-level `if` statements mixed with config blocks.)
    if (params.flair_extra_args?.toString() in ['true', 'false']) {
        throw new IllegalArgumentException(
            "params.flair_extra_args resolved to the literal string '${params.flair_extra_args}' — " +
            "this means --flair_extra_args was followed by a bare flag (e.g. '--flair_extra_args --stringent') " +
            "and Nextflow's CLI parser swallowed it as a separate boolean param instead of this value. " +
            "Use '--flair_extra_args=--stringent' (with '=') instead."
        )
    }

    // ========================================================================================
    // PREPARING INPUT CHANNELS
    // ========================================================================================
    channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->

            def meta = [
                id          : row.sample_id,
                group_id    : row.group_id,
                read_type   : row.read_type,
                library_type: row.library_type,
                replicate   : row.replicate,
                strandedness: row.strandedness,
                single_end  : row.read_type == 'long'
            ]

            def reads = []

            if (row.fastq_1) {
                reads = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
            }

            def bam     = row.bam     ? file(row.bam)     : []
            def bai     = row.bai     ? file(row.bai)     : []
            def star_sj = row.star_sj ? file(row.star_sj) : []
            def gtf     = row.gtf     ? file(row.gtf)     : []
            def bed     = row.bed     ? file(row.bed)     : []
            def fa      = row.fa      ? file(row.fa)      : []
            def map     = row.map     ? file(row.map)     : []

            if (meta.read_type == 'long' && params.skip_aligning_minimap2 && (!bam || !bai) && !params.skip_transcriptome_flair) {
                error "Missing bam/bai for long-read sample '${meta.id}' while --skip_aligning_minimap2 is true"
            }

            if (meta.read_type == 'short' && params.skip_aligning_star && !star_sj) {
                error "Missing star_sj for short-read sample '${meta.id}' while --skip_aligning_star is true"
            }

            if (meta.read_type == 'long' && !reads) {
                error "Missing fastq_1 for long-read sample '${meta.id}'"
            }

            if (meta.read_type == 'short' && !params.skip_aligning_star && reads.size() < 2) {
                error "Missing fastq_1/fastq_2 for short-read sample '${meta.id}'"
            }

            tuple(meta, reads, bam, bai, star_sj, gtf, bed, fa, map)
        }
        .multiMap { meta, reads, bam, bai, star_sj, gtf, bed, fa, map ->
            long_fastq  : tuple(meta, reads)
            long_bam    : tuple(meta, bam, bai)
            short_fastq : tuple(meta, reads)
            short_sj    : tuple(meta, star_sj)
            prebuilt    : tuple(meta, gtf, bed, fa, map)
        }
        .set { ch_input }


    // Split into SHORT and LONG right away
    long_reads_ch  = ch_input.long_fastq.filter  { meta, reads -> meta.read_type == 'long' && reads }
    short_reads_ch = ch_input.short_fastq.filter { meta, reads -> meta.read_type == 'short' && reads }


    genome_fa_ch = channel.fromPath(params.genome).map { fasta ->
            tuple([id: "ref"], fasta)
        }
        .first()

    annotation_ch = channel.fromPath(params.annotation).map { gtf ->
            tuple([id: "ref"], gtf)
        }
        .first()


    
    




    // ========================================================================================
    // PREPROCESS_LONG_READS
    // ========================================================================================
    if (!params.skip_long_preprocessing) {
        preprocess_long_out = PREPROCESS_LONG_READS( long_reads_ch )
        long_reads_for_mapping = preprocess_long_out.filtered_reads
        
    } else {
        long_reads_for_mapping = long_reads_ch.map { meta, reads ->
            tuple(meta, reads[0])
        }
    }





    // ========================================================================================
    // PREPROCESS_SHORT_READS
    // ========================================================================================
    if (!params.skip_short_preprocessing) {
        preprocess_short_out = PREPROCESS_SHORT_READS( short_reads_ch )
        short_reads_for_mapping = preprocess_short_out.short_reads
    } else {
        short_reads_for_mapping = short_reads_ch.map { meta, reads ->
            tuple(meta, reads)
        }    
    }


    // ========================================================================================
    // SAMTOOLS FAIDX
    // ========================================================================================
    SAMTOOLS_FAIDX( genome_fa_ch )




    // ========================================================================================
    // MINIMAP2
    // ========================================================================================

    if (!params.skip_aligning_minimap2) {
     
        minimap2_out = MINIMAP2(
            long_reads_for_mapping,
            genome_fa_ch
        )

        bam_with_bai = minimap2_out.bam
            .join(minimap2_out.bai)
            .map { meta, bam, bai ->
                tuple(meta, bam, bai)
            }

    } else {

        bam_with_bai = ch_input.long_bam.filter { meta, bam, bai -> bam && bai }

    }



    // ========================================================================================
    // SHORT-READ SPLICE JUNCTIONS, OPTIONAL PER tissue_type
    // ========================================================================================

    if (!params.skip_aligning_star) {

        STAR_JUNCTIONS(
            short_reads_for_mapping,
            genome_fa_ch,
            annotation_ch
        )

        star_sj_by_group = STAR_JUNCTIONS.out.short_junctions
            .map { meta, splice_junctions ->
                tuple(meta.group_id, splice_junctions)
            }
            .groupTuple(by: 0)

    } else {

        star_sj_by_group = ch_input.short_sj
            .filter { meta, sj -> meta.read_type == 'short' && sj }
            .map { meta, sj ->
                tuple(meta.group_id, sj)
            }
            .groupTuple(by: 0)
        
    }


    



    // ========================================================================================
    // INTRONPROSPECTOR JUNCTIONS + FLAIR TRANSCRIPTOME
    // ========================================================================================
    if (!params.skip_transcriptome_flair) {

        int_out = INTRONPROSPECTOR_JUNCTIONS(
            bam_with_bai,
            SAMTOOLS_FAIDX.out.fa,
            SAMTOOLS_FAIDX.out.fai,
        )

        flair_out = FLAIR_TRANSCRIPTOME( 
            bam_with_bai,
            int_out.junc_bed,
            star_sj_by_group,
            SAMTOOLS_FAIDX.out.fa,
            annotation_ch
        )
    } else {

        def _flair_raw = ch_input.prebuilt
            // tuple(meta, gtf, bed, fa, map)

        flair_out = [
            flair_gtf: _flair_raw.map { meta, gtf, bed, fa, map -> tuple(meta, gtf) },
            flair_bed: _flair_raw.map { meta, gtf, bed, fa, map -> tuple(meta, bed) },
            flair_fa:  _flair_raw.map { meta, gtf, bed, fa, map -> tuple(meta, fa)  },
            flair_map: _flair_raw.map { meta, gtf, bed, fa, map -> tuple(meta, map) }
        ]
    }




    // ========================================================================================
    // FLAIR COMBINE QUANTIFY
    // ========================================================================================
    flair_cq_out = FLAIR_COMBINE_QUANTIFY( 
        flair_out.flair_bed,
        flair_out.flair_fa, 
        flair_out.flair_map, 
        flair_out.flair_gtf, 
        long_reads_for_mapping
    )




    // ========================================================================================
    // SQANTI
    // ========================================================================================
    sqanti_out = SQANTI(
        flair_cq_out.combined_gtf,
        SAMTOOLS_FAIDX.out.fa,
        annotation_ch,
        flair_cq_out.quantify_counts
    )




    // ========================================================================================
    // FINDING OVERLAPS
    // ========================================================================================
    AGAT_OVERLAP(
        sqanti_out.sqanti_gtf_rescued,
        sqanti_out.sqanti_classification
    )

}
