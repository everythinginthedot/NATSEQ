include { STAR_ALIGN          } from '../../../modules/nf-core/star/align/main'  
include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate/main'


workflow STAR_JUNCTIONS {
    take:
    short_reads_ch   // tuple val(meta), [path(fastq_1), path(fastq_2)]
    genome_fa
    annotation


    main:

    annotation = annotation
        .map{ gtf -> tuple([id: "ref"], gtf) }


    // STAR index — built once, shared across all samples
    if (params.star_index) {

        star_index = channel.fromPath(params.star_index)
            .map { star ->
                tuple([id: "star_index"], star)
            }
            .first()

    } else {

        STAR_GENOMEGENERATE(
            genome_fa,
            annotation
        )

        star_index = STAR_GENOMEGENERATE.out.index
    }



    // Align short reads
    STAR_ALIGN(
        short_reads_ch,
        star_index,
        annotation,
        params.star_ignore_sjdbgtf
    )



    emit:
    short_junctions = STAR_ALIGN.out.spl_junc_tab  
    // tuple(meta, SJ.out.tab); downstream groups junctions by meta.group_id
  
}