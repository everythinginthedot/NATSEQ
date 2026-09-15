process FLAIR {

    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "brookslab/flair:3.0.0"

    input:
    tuple val(meta), path(bam), path(bai), path(junctions), path(junctions_short)   // short read SJ.out.tab 
    tuple val(meta2), path(genome_fa)
    tuple val(meta3), path(gtf)
    

    output:
    tuple val(meta), path("*.log"),        emit: log
    tuple val(meta), path("*.counts.txt"), emit: counts
    tuple val(meta), path("*.map.txt"),    emit: map
    tuple val(meta), path("*.bed"),        emit: bed
    tuple val(meta), path("*.fa"),         emit: fa
    tuple val(meta), path("*.gtf"),        emit: gtf


    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def junc_arg = ""

    if ( junctions_short.toString() != "[]" && junctions_short ) {
        junc_arg = "--junction_tab ${junctions_short}"
    } 
    // 2. No Illumina — fall back to the long-read BED
    else if ( junctions.toString() != "[]" && junctions ) {
        junc_arg = "--junction_bed ${junctions}"
    }

    """   
    flair transcriptome \
        -g ${genome_fa} \
        -f ${gtf} \
        ${junc_arg} \
        -b ${bam} \
        -o ${prefix}.flair \
        ${args} \
        > ${prefix}.flair.log 2>&1
    """
}

//         

// flair transcriptome         -g data/ref/GRCh38.primary_assembly.genome.fa         -f data/gtf/gencode.v49.annotation.gtf         --junction_bed results/intron_junctions/SRR30264742.IPjunctions.bed         --junction_support 2         -b results/bam/SRR30264740.bam         -o SRR30264740.flair         > SRR30264740.flair.log 2>&1
// flair transcriptome         -g data/ref/GRCh38.primary_assembly.genome.fa         -f data/gtf/gencode.v49.annotation.gtf         --junction_bed results/intron_junctions/SRR30264742.IPjunctions.bed         --junction_support 2         -b results/bam/SRR30264740/SRR30264740.bam         -o SRR30264740.flair         > SRR30264740.flair.log 2>&1