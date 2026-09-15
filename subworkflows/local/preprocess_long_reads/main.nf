/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { NANOPLOT as NANOPLOT_BEFORE   } from '../../../modules/nf-core/nanoplot/main'     
include { NANOPLOT as NANOPLOT_AFTER    } from '../../../modules/nf-core/nanoplot/main'    
include { PORECHOP_PORECHOP as PORECHOP } from '../../../modules/nf-core/porechop/porechop/main'       
include { PORECHOP_ABI                  } from '../../../modules/nf-core/porechop/abi/main'                                                                                            
include { FILTLONG                      } from '../../../modules/nf-core/filtlong/main'        
include { NANOFILT                      } from '../../../modules/nf-core/nanofilt/main'    
include { CHOPPER                       } from '../../../modules/nf-core/chopper/main'    
include { SEQKIT_SPLIT2                 } from '../../../modules/nf-core/seqkit/split2/main'
include { FASTQ_CONCAT                  } from '../../../modules/local/cat/main'
include { PYCHOPPER                     } from '../../../modules/nf-core/pychopper/main'                                                                               


workflow PREPROCESS_LONG_READS {
    take:
    long_reads_ch


    main:

    // 1. Initial QC of the raw (whole) input file
    NANOPLOT_BEFORE( long_reads_ch )



    // ========================================================================================
    // STEP 1: DYNAMIC CHUNK-SIZE CALCULATION AND SPLITTING (GLOBAL, ON INPUT)
    // ========================================================================================
    long max_monolith_bytes = MemoryUnit.of(params.max_monolith_size).toBytes()
    long target_chunk_bytes  = MemoryUnit.of(params.target_chunk_size).toBytes()

    ch_evaluated_reads = long_reads_ch.map { meta, fastq ->
        def file_size = fastq[0].size()
        def chunks = 1
        if (file_size > max_monolith_bytes) {
            chunks = Math.ceil(file_size / target_chunk_bytes).toInteger()
        }
        return tuple(meta + [split_chunks: chunks], fastq)
    } // add split_chunks to meta

    // Split the stream: large files go to SeqKit, small ones pass through as-is
    ch_evaluated_reads.branch { meta, fastq ->
        split:   meta.split_chunks > 1
        monolith: true
    }.set { ch_split_logic }

    // Split the large files
    SEQKIT_SPLIT2( ch_split_logic.split )

    // Flatten the resulting chunks into one stream, giving each its own meta.id
    ch_prepared_chunks = SEQKIT_SPLIT2.out.reads
        .transpose()
        .map { meta, chunk ->
            def m = (chunk.name =~ /part[_.]?(\d+)/)
            def part = m ? m[0][1] : '0'
            def chunk_meta = meta + [
                id: "${meta.id}_part${part}",
                original_meta: meta,
                part: part
            ]
            tuple(chunk_meta, chunk)
        }

    // Merge monoliths and split chunks back into a single channel — every
    // downstream tool just sees "a stream of similarly-sized FASTQ files".
    ch_processing_stream = ch_prepared_chunks.mix( ch_split_logic.monolith )






    // ========================================================================================
    // STEP 2: ADAPTER TRIMMING (STREAMING)
    // ========================================================================================
    // Tool is picked per-sample from meta.library_type, not one global
    // params.adaptertrimming_tool for the whole run — lets a samplesheet mix
    // dRNA/cDNA/direct-cDNA safely. Pychopper also reorients reads to the correct
    // strand by primer identification, which porechop does not do.
    ch_tool_routed = ch_processing_stream.map { meta, fastq ->
        def tool
        if (meta.library_type == 'dRNA') {
            tool = params.adaptertrimming_tool_drna
        } else if (meta.library_type == 'cDNA') {
            tool = params.adaptertrimming_tool_cdna
        } else if (meta.library_type == 'direct-cDNA') {
            tool = params.adaptertrimming_tool_directcdna
        } else {
            tool = params.adaptertrimming_tool
        }
        tuple(meta + [adaptertrimming_tool_resolved: tool], fastq)
    }

    ch_tool_routed.branch { meta, fastq ->
        porechop:     meta.adaptertrimming_tool_resolved == 'porechop'
        porechop_abi: meta.adaptertrimming_tool_resolved == 'porechop_abi'
        pychopper:    meta.adaptertrimming_tool_resolved == 'pychopper'
        other:        true
    }.set { ch_by_tool }

    // Fail loudly on an unrecognized value — same as the old global error(), just
    // checked per-sample now instead of for the whole run at once.
    ch_by_tool.other.subscribe { meta, fastq ->
        error "Unsupported adaptertrimming_tool '${meta.adaptertrimming_tool_resolved}' for sample '${meta.id}' (library_type='${meta.library_type}')"
    }

    PORECHOP( ch_by_tool.porechop.map { meta, fastq -> tuple(meta, fastq) } )
    PORECHOP_ABI( ch_by_tool.porechop_abi.map { meta, fastq -> tuple(meta, fastq) }, params.custom_adapters_porechop ?: [] )
    PYCHOPPER( ch_by_tool.pychopper.map { meta, fastq -> tuple(meta, fastq) } )

    // PYCHOPPER has no separate log channel to mix in (same as before).
    trimmed_reads   = PORECHOP.out.reads.mix( PORECHOP_ABI.out.reads ).mix( PYCHOPPER.out.fastq )
    trimmed_log     = PORECHOP.out.log.mix( PORECHOP_ABI.out.log )
    trimmed_version = PORECHOP.out.versions.mix( PORECHOP_ABI.out.versions ).mix( PYCHOPPER.out.versions )

    



    // ========================================================================================
    // STEP 3: QUALITY FILTERING (SAME STREAM)
    // ========================================================================================

    if (params.filtering_tool == 'chopper') {
        CHOPPER( trimmed_reads, [] )
        filtered_stream  = CHOPPER.out.fastq
        filtered_log     = ''
        filtered_version = CHOPPER.out.versions_chopper
    } else if (params.filtering_tool == 'filtlong') {
        FILTLONG( trimmed_reads.map { meta, lr -> tuple(meta, [], lr) } )
        filtered_stream  = FILTLONG.out.reads
        filtered_log     = FILTLONG.out.log
        filtered_version = FILTLONG.out.versions_filtlong
    } else if (params.filtering_tool == 'nanofilt') {
        NANOFILT( trimmed_reads, [] )
        filtered_stream  = NANOFILT.out.filtreads
        filtered_log     = NANOFILT.out.log_file
        filtered_version = NANOFILT.out.versions
    } else {
        error "Unsupported filtering_tool: ${params.filtering_tool}"
    }




    // ========================================================================================
    // STEP 4: REASSEMBLE CHUNKS BACK INTO MONOLITHS (GLOBAL, ON OUTPUT)
    // ========================================================================================

    // Split the stream again: former chunks get grouped for concatenation, former
    // monoliths pass through unchanged.
    filtered_stream.branch { meta, fastq ->
        to_concat: meta.split_chunks > 1
        as_is:     true
    }.set { ch_post_logic }

    // Group chunks by their original meta (restoring original_meta). groupKey(id,
    // expected_chunk_count) instead of a bare id — otherwise groupTuple(by: 0)
    // without size can't tell a sample's group is complete until the WHOLE channel
    // closes (i.e. until every sample's chunks in the run have been filtered) — same
    // bug class already fixed in the FLAIR/SQANTI subworkflows. The expected chunk
    // count is already known (meta.original_meta.split_chunks, computed in Step 1),
    // so groupKey can be used exactly, without a separate size:.
    ch_grouped_chunks = ch_post_logic.to_concat
        .map { meta, reads ->
            tuple(groupKey(meta.original_meta.id, meta.original_meta.split_chunks), meta.original_meta, reads)
        }
        .groupTuple(by: 0)
        .map { original_id, metas, reads ->
            tuple(metas[0], reads.flatten())
        }

    // Concatenate grouped chunks into the final .fastq.gz files
    FASTQ_CONCAT( ch_grouped_chunks )

    // Merge concatenated files back with the ones that were always monolithic
    filtered_reads = FASTQ_CONCAT.out.reads.mix( ch_post_logic.as_is )

    // 2. Final QC of the reassembled file
    NANOPLOT_AFTER ( filtered_reads )




    emit:

    be_nanoplot_png     = NANOPLOT_BEFORE.out.png
    be_nanoplot_html    = NANOPLOT_BEFORE.out.html
    be_nanoplot_txt     = NANOPLOT_BEFORE.out.txt
    be_nanoplot_log     = NANOPLOT_BEFORE.out.log
    be_nanoplot_version = NANOPLOT_BEFORE.out.versions

    trimmed_reads
    trimmed_log
    trimmed_version

    filtered_reads
    filtered_log
    filtered_version

    af_nanoplot_png     = NANOPLOT_AFTER.out.png
    af_nanoplot_html    = NANOPLOT_AFTER.out.html
    af_nanoplot_txt     = NANOPLOT_AFTER.out.txt
    af_nanoplot_log     = NANOPLOT_AFTER.out.log
    af_nanoplot_version = NANOPLOT_AFTER.out.versions
}
