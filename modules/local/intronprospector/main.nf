process INTRONPROSPECTOR {

    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "juanj24/intronprospector:latest"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(genome_fa) 
    tuple val(meta3), path(genome_fai)

    
    

    output:
    tuple val(meta), path("*.bed"), emit: junctions

    script:
    """
    export LD_LIBRARY_PATH="\${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH:-}"

    intronProspector -C 0.0 \
        ${bam} \
        --intron-bed6=${bam.baseName}.IPjunctions.bed \
        --genome-fasta=${genome_fa}
    """
}