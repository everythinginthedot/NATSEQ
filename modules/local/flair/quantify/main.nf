process FLAIR_QUANTIFY {
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "brookslab/flair:3.0.0"

    input:
    tuple val(meta), val(samples), path(fastqs), path(comb_fa)

    output:
    tuple val(meta), path("*.counts.tsv"), emit: counts

    script:
    def args = task.ext.args ?: ""

    def manifest_lines = [samples, fastqs]
        .transpose()
        .collect { sample_meta, fastq ->
            def clean_id    = sample_meta.id.replaceAll('_', '-') // FLAIR's manifest format uses '_' as an internal delimiter
            def cond        = sample_meta.condition ?: 'condition1'
            def batch       = sample_meta.batch     ?: 'batch1'
            def clean_cond  = cond.replaceAll('_', '-')
            def clean_batch = batch.replaceAll('_', '-')


            "${clean_id}\t${clean_cond}\t${clean_batch}\t${fastq}"
        }
        .join('\n')

    """
    cat > quantify_manifest.tsv << 'EOF'
${manifest_lines}
EOF

    flair quantify \
        ${args} \
        --reads_manifest quantify_manifest.tsv \
        --isoforms ${comb_fa} \
        --output ${meta.id}.counts
    """
}