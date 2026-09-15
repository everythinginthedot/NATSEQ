process SQANTI3_QC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    // gtf + flair_counts MUST travel together in one tuple, keyed by the same
    // meta — they come from two independent upstream channels (combined_gtf vs
    // quantify_counts) that don't necessarily emit same-sample items in the same
    // relative order. Two separate positional inputs would pair them by arrival
    // order, not by meta.id, and silently mismatch samples (this happened once —
    // seed3_0.05's GTF got paired with seed3_0.50's counts).
    tuple val(meta), path(gtf), path(flair_counts)
    tuple val(meta2), path(ref_fasta)
    tuple val(meta3), path(ref_gtf)


    output:
    tuple val(meta), path("${prefix}_QC_output/"),                     emit: results     // whole output dir
    tuple val(meta), path("${prefix}_QC_output/*.txt"),                emit: txt,            optional: true
    tuple val(meta), path("${prefix}_QC_output/${prefix}_classification.txt"),  emit: classification, optional: true
    tuple val(meta), path("${prefix}_QC_output/*.pdf"),                emit: pdf,            optional: true
    tuple val(meta), path("${prefix}_QC_output/*.tsv"),                emit: tsv,            optional: true
    tuple val(meta), path("${prefix}_QC_output/${prefix}_corrected.fasta"),              emit: fasta,          optional: true
    tuple val(meta), path("${prefix}_QC_output/${prefix}_corrected.gtf"),      emit: gtf,            optional: true

    path "versions.yml",                                               emit: versions


    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export PATH="${projectDir}/bin/SQANTI3:\$PATH"

    # FLAIR_QUANTIFY's counts.tsv needs two fixes before SQANTI3's FLcount_parser
    # will read it correctly:
    # 1. Header/delimiter: it ships as "ids<TAB>sample1[<TAB>sample2...]", but the
    #    parser wants either a single-sample header "pbid<TAB>count_fl", or a
    #    multi-sample header starting with "superPBID" (tab-separated) or "id"
    #    (but THAT variant is comma-separated, not tab — a literal "id"-prefixed
    #    tab file gets misparsed).
    # 2. ID format: FLAIR's row keys are "<transcript_id>_<gene_id>" (e.g.
    #    "ENST00000037502.11_ENSG00000034971.18", or "..._chr1:12345000" for
    #    unannotated loci), while the --isoforms GTF's transcript_id is plain.
    #    Verified on real data: 0/6910 exact matches before stripping the gene
    #    suffix, all matched after — without this every isoform silently gets
    #    FL=0 regardless of real read support (found 2026-08-24, see project memory).
    #    Strip using the GTF's own transcript_id->gene_id mapping rather than a
    #    regex guess at what a gene_id can look like.
    python3 - "${gtf}" "${flair_counts}" fl_counts_for_sqanti.tsv <<'PYEOF'
import re, sys

gtf_path, counts_path, out_path = sys.argv[1:4]

tx2gene = {}
with open(gtf_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\\n').split('\\t')
        if len(fields) < 9 or fields[2] != 'transcript':
            continue
        tx_m = re.search(r'transcript_id "([^"]+)"', fields[8])
        gene_m = re.search(r'gene_id "([^"]+)"', fields[8])
        if tx_m and gene_m:
            tx2gene[tx_m.group(1)] = gene_m.group(1)
compound_to_plain = {f"{tx}_{gene}": tx for tx, gene in tx2gene.items()}

with open(counts_path) as fin, open(out_path, 'w') as fout:
    header = fin.readline().rstrip('\\n').split('\\t')
    n_samples = len(header) - 1
    fout.write(('pbid\\tcount_fl\\n' if n_samples == 1 else 'superPBID\\t' + '\\t'.join(header[1:]) + '\\n'))
    for line in fin:
        parts = line.rstrip('\\n').split('\\t')
        parts[0] = compound_to_plain.get(parts[0], parts[0])
        fout.write('\\t'.join(parts) + '\\n')
PYEOF

    sqanti3_qc.py \
        ${args} \
        --isoforms ${gtf} \
        --refGTF ${ref_gtf} \
        --refFasta ${ref_fasta} \
        --fl_count fl_counts_for_sqanti.tsv \
        --dir ${prefix}_QC_output \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_qc.py --version 2>&1 | sed 's/SQANTI3 version //')
    END_VERSIONS
    """
}