process FLAIR_COMBINE {
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "brookslab/flair:3.0.0"

    input:
    tuple val(meta), val(samples),
          path(beds,  stageAs: 'input/*'),
          path(fas,   stageAs: 'input/*'),
          path(maps,  stageAs: 'input/*')

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.fa"),  emit: fa
    tuple val(meta), path("*.gtf"), emit: gtf, optional: true
    tuple val(meta), path("*.counts.txt"), emit: counts, optional: true
    tuple val(meta), path("*.map.txt"),    emit: map
    tuple val(meta), path("*.log"),        emit: log, optional: true
    
    script:
    def args   = task.ext.args ?: ""
    // Build the manifest right here — the files are already staged in input/
    def manifest_lines = [samples, beds, fas, maps]
        .transpose()
        .collect { sample_meta, bed, fa, map ->
            "${sample_meta.id}\tisoform\t${bed}\t${fa}\t${map}"
        }
        .join('\n')

    """
    # Manifest with local paths (files are already in the work dir)
    cat > manifest.tsv << 'EOF'
${manifest_lines}
EOF

    flair combine \\
        ${args} \\
        --manifest manifest.tsv \\
        -o ${meta.id}
    """
}