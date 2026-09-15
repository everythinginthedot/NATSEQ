/*
========================================================================================
    IMPORT NF-CORE MODULES
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index/main'  
include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'  
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_MINIMAP2 } from '../../../modules/nf-core/seqkit/split2/main'
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'


workflow MINIMAP2 {

    take:
    reads_ch
    genome_fa_ch

    main:


    // ========================================================================================
    // STEP 1: index the genome
    // ========================================================================================
    MINIMAP2_INDEX(genome_fa_ch)



    // ========================================================================================
    // STEP 2: split files into chunks
    // ========================================================================================
    long max_monolith_bytes = MemoryUnit.of(params.minimap2_max_monolith_size).toBytes()
    long target_chunk_bytes = MemoryUnit.of(params.minimap2_target_chunk_size).toBytes()

    ch_evaluated_reads = reads_ch.map { meta, fastq ->
        def file_size = fastq.size()
        def chunks = 1

        if (file_size > max_monolith_bytes) {
            chunks = Math.ceil(file_size / target_chunk_bytes).toInteger()
        }

        tuple(meta + [minimap2_split_chunks: chunks], fastq)
    }

    ch_evaluated_reads.branch { meta, fastq ->
        split:    meta.minimap2_split_chunks > 1
        monolith: true
    }.set { ch_split_logic }


    SEQKIT_SPLIT2_MINIMAP2(ch_split_logic.split)



    ch_split_reads = SEQKIT_SPLIT2_MINIMAP2.out.reads
        .transpose()
        .map { meta, chunk ->
            def m = (chunk.name =~ /part[_.]?(\d+)/)
            def part = m ? m[0][1] : '0'

            def chunk_meta = meta + [
                id           : "${meta.id}_part${part}",
                original_meta: meta,
                part         : part
            ]

            tuple(chunk_meta, chunk)
        }


    ch_reads_for_align = ch_split_reads.mix(ch_split_logic.monolith)




    // ========================================================================================
    // STEP 3: align to the genome
    // ========================================================================================
    MINIMAP2_ALIGN(
        ch_reads_for_align,
        MINIMAP2_INDEX.out.index,
        true,
        '',
        params.cigar_paf_format,
        params.cigar_bam
    )

    


    // ========================================================================================
    // STEP 4: reassemble chunks back into files
    // ========================================================================================
    MINIMAP2_ALIGN.out.bam
    .branch { meta, bam ->
        chunk:    meta.original_meta != null
        monolith: true
    }
    .set { ch_aligned_bam }


    // groupKey(id, expected_chunk_count) instead of a bare id — otherwise
    // groupTuple(by: 0) without size can't tell a sample's group is complete until
    // the WHOLE channel closes (i.e. until every sample's chunks in the run have
    // finished aligning) — same bug class already fixed in FLAIR/SQANTI and in
    // preprocess_long_reads. The expected chunk count is already known
    // (meta.original_meta.minimap2_split_chunks, computed in Step 2).
    ch_bam_for_merge = ch_aligned_bam.chunk
        .map { meta, bam ->
            def sample_meta = meta.original_meta ?: meta
            tuple(groupKey(sample_meta.id, sample_meta.minimap2_split_chunks), sample_meta, bam)
        }
        .groupTuple(by: 0)
        .map { sample_id, metas, bams ->
            tuple(metas[0], bams, [])
        }


    ch_empty_ref_for_merge = channel.value(
        tuple([id: 'no_ref'], [], [], [])
    )


    SAMTOOLS_MERGE(
        ch_bam_for_merge,
        ch_empty_ref_for_merge
    )


    ch_final_bam = SAMTOOLS_MERGE.out.bam.mix(
        ch_aligned_bam.monolith
    )


    // ========================================================================================
    // STEP 5: index the assembled bam files
    // ========================================================================================
    SAMTOOLS_INDEX(ch_final_bam)



    emit:
    index_ind      = MINIMAP2_INDEX.out.index
    versions_index = MINIMAP2_INDEX.out.versions_minimap2

    versions_align = MINIMAP2_ALIGN.out.versions_minimap2

    bam            = ch_final_bam
    bai            = SAMTOOLS_INDEX.out.index
}