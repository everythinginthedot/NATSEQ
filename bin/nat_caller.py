#!/usr/bin/env python3
"""
NAT scanner for reconstructed transcript annotations.

The primary analysis scans one query annotation against itself and reports
opposite-strand transcript pairs that are both present in the reconstructed
GTF/GFF. A reference annotation is used as an annotation/confidence layer.

The secondary analysis reports weaker reference-anchored antisense hints:
single query transcripts that overlap opposite-strand reference transcripts.
"""

import argparse
import csv
import json
import os
import re
import sqlite3
import sys
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass

try:
    import gffutils
except ModuleNotFoundError:
    gffutils = None


__version__ = "1.1.0-l0"

TRANSCRIPT_TYPES = ("transcript", "mRNA")
EXON_TYPES = ("exon",)

# --- L0/curation additions: biotype + tag vocabulary ----------------------
# GENCODE writes gene_type/transcript_type; Ensembl writes *_biotype.
GENE_TYPE_KEYS = ("gene_type", "gene_biotype", "biotype")
TX_TYPE_KEYS = ("transcript_type", "transcript_biotype")
TAG_KEYS = ("tag",)

# Curation classes follow the convention used in the golden dataset:
#   p = protein coding, n = non-coding RNA, o = other (pseudogene, IG/TR),
#   u = unknown (attribute absent, e.g. de novo FLAIR transcripts)
PROTEIN_CODING_TYPES = {
    "protein_coding",
    "protein_coding_LoF",
    "protein_coding_CDS_not_defined",
}
# Genes that are typically processed CO-TRANSCRIPTIONALLY from a host
# transcript's intron rather than transcribed independently from their own
# promoter (Brown & Ares 1993 for snoRNA; Baskerville & Bartel 2005 for
# intronic miRNA; scaRNA shares the snoRNA processing pathway). Such genes
# are not expected to appear as an independent long-read transcript, so
# antisense pairs involving them measure coincidence with the host's
# coordinates, not the pipeline's actual detection ability.
# NOTE: this is a biotype-level proxy, not a per-locus verified check - some
# loci in these classes do have independent promoters (e.g. intergenic
# miRNAs). Kept deliberately narrow: snRNA is EXCLUDED because the major
# spliceosomal snRNAs (U1/U2/U4/U5) are independently transcribed.
COTRANSCRIPTIONAL_SMALL_RNA_TYPES = {"miRNA", "snoRNA", "scaRNA"}

NONCODING_TYPES = {
    "lncRNA", "lincRNA", "antisense", "antisense_RNA", "processed_transcript",
    "sense_intronic", "sense_overlapping", "bidirectional_promoter_lncRNA",
    "macro_lncRNA", "3prime_overlapping_ncRNA", "non_coding", "miRNA",
    "snRNA", "snoRNA", "scaRNA", "misc_RNA", "rRNA", "Mt_tRNA", "Mt_rRNA",
    "ribozyme", "sRNA", "scRNA", "vault_RNA", "TEC",
}

# Overlap-length strata used by the benchmark (bp). Closed-open on the left,
# closed on the right: [1,100) / [100,500] / (500, inf)
OVERLAP_STRATUM_BREAKS = (100, 500)


@dataclass(frozen=True)
class TxRecord:
    transcript_id: str
    gene_id: str
    gene_name: str
    seqid: str
    start: int
    end: int
    strand: str
    exons: tuple
    attrs: dict
    # --- L0/curation additions (defaults keep old constructors working) ---
    gene_type: str = ""
    transcript_type: str = ""
    tags: tuple = ()

    @property
    def span_len(self):
        return self.end - self.start + 1

    @property
    def exonic_len(self):
        return interval_length(self.exons)

    @property
    def n_exons(self):
        return len(self.exons)

    @property
    def biotype_class(self):
        return classify_biotype(self.gene_type)

    @property
    def is_mane(self):
        return any(str(t).startswith("MANE_Select") for t in self.tags)

    @property
    def is_canonical(self):
        return "Ensembl_canonical" in self.tags

    @property
    def is_basic(self):
        return "basic" in self.tags

    @property
    def is_appris_principal(self):
        return any(str(t).startswith("appris_principal") for t in self.tags)

    @property
    def tsl(self):
        """transcript_support_level, first token only ('1 (assigned...)' -> '1')."""
        raw = self.attrs.get("transcript_support_level", "")
        return str(raw).split()[0] if raw else ""

    @property
    def is_tagged(self):
        """Carries at least one 'representative isoform' tag (see DR-02)."""
        return self.is_mane or self.is_canonical or self.is_appris_principal


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scan reconstructed GTF/GFF annotations for candidate natural antisense transcript pairs."
    )
    parser.add_argument("--version", action="version", version=f"nat_caller.py {__version__}")
    parser.add_argument("--query", required=True, help="Query GTF/GFF, e.g. FLAIR/SQANTI transcript annotation.")
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference GTF/GFF annotation, e.g. GENCODE. Used for annotation and weaker reference hints.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output TSV for reconstructed query-vs-query NAT pairs.")
    parser.add_argument(
        "--reference-hints-output",
        default=None,
        help="Output TSV for weaker query-vs-reference antisense hints. Defaults to '<output>.reference_hints.tsv'.",
    )
    parser.add_argument(
        "--sqanti-classification",
        default=None,
        help="Optional SQANTI classification TSV/TXT. If supplied, query records are annotated with SQANTI columns.",
    )
    parser.add_argument(
        "--min-first-exon-overlap-bp",
        type=int,
        default=0,
        help="Optional minimum overlap between the first exons of both transcripts. Default 0 keeps span-only candidates.",
    )
    parser.add_argument("--min-span-overlap-bp", type=int, default=1, help="Minimum transcript-span overlap in bp.")
    parser.add_argument(
        "--min-query-exonic-fraction",
        type=float,
        default=0.0,
        help="Minimum fraction of first transcript/query exonic bases overlapped.",
    )
    parser.add_argument(
        "--min-reference-exonic-fraction",
        type=float,
        default=0.0,
        help="Minimum fraction of second transcript/reference exonic bases overlapped.",
    )
    parser.add_argument(
        "--legacy-anchoring",
        action="store_true",
        help="Reproduce pre-1.1.0 reference anchoring (see best_reference_hit docstring).",
    )
    parser.add_argument(
        "--keep-db",
        action="store_true",
        help="Keep temporary gffutils SQLite databases next to the output.",
    )
    parser.add_argument(
        "--tmpdir",
        default=None,
        help="Directory for temporary gffutils databases. Defaults to system temp.",
    )
    return parser.parse_args()


def log(message):
    print(f"[nat_caller] {message}", file=sys.stderr, flush=True)


def get_attr(feature, keys, default=""):
    if isinstance(keys, str):
        keys = (keys,)
    for key in keys:
        values = feature.attributes.get(key)
        if values:
            return values[0]
    return default


def attr_dict(feature):
    return {key: values[0] if values else "" for key, values in feature.attributes.items()}


def attr_list(feature, keys):
    """Return ALL values of the first present key.

    Needed because GENCODE repeats the 'tag' attribute (basic, MANE_Select,
    Ensembl_canonical, appris_principal_1, ...) and attr_dict() keeps only the
    first value, silently discarding the rest.
    """
    if isinstance(keys, str):
        keys = (keys,)
    for key in keys:
        values = feature.attributes.get(key)
        if values:
            return tuple(values)
    return tuple()


def classify_biotype(gene_type):
    """Map a GENCODE/Ensembl gene_type onto the curation class p/n/o/u."""
    if not gene_type:
        return "u"
    if gene_type in PROTEIN_CODING_TYPES:
        return "p"
    if gene_type in NONCODING_TYPES:
        return "n"
    return "o"


def biotype_pair_class(class_a, class_b):
    """Order-independent pair label: p-p, n-p, n-n, o-p, ..."""
    return "-".join(sorted((class_a or "u", class_b or "u")))


def overlap_length_stratum(bp):
    if bp <= 0:
        return "none"
    if bp < OVERLAP_STRATUM_BREAKS[0]:
        return "lt100"
    if bp <= OVERLAP_STRATUM_BREAKS[1]:
        return "100_500"
    return "gt500"


def igv_locus(tx1, tx2, pad=500):
    """IGV-pasteable locus covering both partners plus flanking context."""
    start = max(1, min(tx1.start, tx2.start) - pad)
    end = max(tx1.end, tx2.end) + pad
    return f"{tx1.seqid}:{start}-{end}"


def sanitize_db_prefix(path):
    name = os.path.basename(path)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def create_feature_db(annotation, output, label, tmpdir=None, keep_db=False):
    if gffutils is None:
        raise SystemExit(
            "ERROR: Python package 'gffutils' is required. Install it in the process environment."
        )

    if keep_db:
        db_path = f"{output}.{label}.db"
    else:
        handle = tempfile.NamedTemporaryFile(
            prefix=f"nat_caller.{label}.",
            suffix=f".{sanitize_db_prefix(annotation)}.db",
            dir=tmpdir,
            delete=False,
        )
        db_path = handle.name
        handle.close()

    log(f"Creating {label} database: {db_path}")
    gffutils.create_db(
        annotation,
        db_path,
        force=True,
        keep_order=True,
        sort_attribute_values=True,
        merge_strategy="merge",
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )
    return db_path, gffutils.FeatureDB(db_path, keep_order=True)


def children_of_type(db, parent, featuretypes):
    out = []
    for featuretype in featuretypes:
        out.extend(db.children(parent, featuretype=featuretype, order_by="start"))
    return out


def feature_id(feature):
    return get_attr(feature, ("transcript_id", "ID", "Name"), feature.id)


def gene_id_for_transcript(tx):
    return get_attr(tx, ("gene_id", "Parent", "gene", "ref_gene_id"), "")


def gene_name_for_transcript(tx, gene_id):
    name = get_attr(tx, ("gene_name", "gene", "ref_gene_name"), "")
    if name and name != gene_id:
        return name
    return ""


def normalize_intervals(intervals):
    if not intervals:
        return tuple()
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        last = merged[-1]
        if start <= last[1] + 1:
            last[1] = max(last[1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def interval_length(intervals):
    return sum(end - start + 1 for start, end in intervals)


def overlap_bp(intervals_a, intervals_b):
    i = 0
    j = 0
    total = 0
    a = normalize_intervals(intervals_a)
    b = normalize_intervals(intervals_b)
    while i < len(a) and j < len(b):
        start = max(a[i][0], b[j][0])
        end = min(a[i][1], b[j][1])
        if start <= end:
            total += end - start + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def first_exon(tx):
    if not tx.exons:
        return tuple()
    if tx.strand == "+":
        return (tx.exons[0],)
    if tx.strand == "-":
        return (tx.exons[-1],)
    return tuple()


def last_exon(tx):
    if not tx.exons:
        return tuple()
    if tx.strand == "+":
        return (tx.exons[-1],)
    if tx.strand == "-":
        return (tx.exons[0],)
    return tuple()


def span_overlap_bp(a, b):
    if a.seqid != b.seqid:
        return 0
    start = max(a.start, b.start)
    end = min(a.end, b.end)
    return max(0, end - start + 1)


def tss(tx):
    return tx.start if tx.strand == "+" else tx.end


def tes(tx):
    return tx.end if tx.strand == "+" else tx.start


def tx_contains(a, b):
    return a.start <= b.start and a.end >= b.end


def classify_orientation(tx1, tx2):
    if tx_contains(tx2, tx1):
        return "embedded_tx1_in_tx2"
    if tx_contains(tx1, tx2):
        return "embedded_tx2_in_tx1"

    if tx1.strand == "+" and tx2.strand == "-":
        return "tail_to_tail" if tx1.start < tx2.start else "head_to_head"
    if tx1.strand == "-" and tx2.strand == "+":
        return "tail_to_tail" if tx2.start < tx1.start else "head_to_head"

    raise ValueError(
        "Cannot classify orientation for pair "
        f"{tx1.transcript_id}({tx1.strand}) and {tx2.transcript_id}({tx2.strand})"
    )


def container_contained(tx1, tx2):
    """Canonical embedded roles, independent of which transcript is 'tx1'.

    classify_orientation() labels embedding relative to argument order, and the
    argument order depends on the sort order of the transcript list. These
    columns are order-free so gene-level aggregation stays deterministic.
    Returns (container_id, contained_id, reciprocal_flag).
    """
    c12 = tx_contains(tx1, tx2)
    c21 = tx_contains(tx2, tx1)
    if c12 and c21:
        return tx1.transcript_id, tx2.transcript_id, "1"
    if c12:
        return tx1.transcript_id, tx2.transcript_id, "0"
    if c21:
        return tx2.transcript_id, tx1.transcript_id, "0"
    return "", "", "0"


def pair_support_class(tx1_anchor, tx2_anchor):
    """Fills the previously-empty pair_support column.

    'same_anchor_gene' is the important one: both partners anchor to the SAME
    reference gene, which is the signature of a strand-assignment artefact
    (sense leakage) rather than a genuine antisense pair.
    """
    gene1 = tx1_anchor["anchor_gene_id"]
    gene2 = tx2_anchor["anchor_gene_id"]
    anchored1 = tx1_anchor["anchor_source"] != "query_gtf"
    anchored2 = tx2_anchor["anchor_source"] != "query_gtf"
    if anchored1 and anchored2:
        if gene1 and gene1 == gene2:
            return "same_anchor_gene"
        return "distinct_anchor_genes"
    if anchored1 or anchored2:
        return "one_anchor_missing"
    return "no_anchor"


def overlap_evidence_type(span_bp, first_exon_bp, last_exon_bp, all_exon_bp):
    evidence = []
    if first_exon_bp > 0:
        evidence.append("first_exon_exon_overlap")
    if last_exon_bp > 0:
        evidence.append("last_exon_exon_overlap")
    if all_exon_bp > 0 and first_exon_bp == 0 and last_exon_bp == 0:
        evidence.append("internal_or_mixed_exon_overlap")
    if evidence:
        return ";".join(evidence)
    if span_bp > 0:
        return "transcript_span_overlap"
    return "none"


def load_transcripts(db, label):
    transcripts = []
    skipped_no_exons = 0

    for featuretype in TRANSCRIPT_TYPES:
        for tx in db.features_of_type(featuretype, order_by=("seqid", "start")):
            tid = feature_id(tx)
            gid = gene_id_for_transcript(tx)
            gname = gene_name_for_transcript(tx, gid or tid)
            gtype = get_attr(tx, GENE_TYPE_KEYS, "")
            ttype = get_attr(tx, TX_TYPE_KEYS, "")
            tags = attr_list(tx, TAG_KEYS)
            exons = children_of_type(db, tx, EXON_TYPES)

            if exons:
                exon_intervals = normalize_intervals((exon.start, exon.end) for exon in exons)
            else:
                exon_intervals = tuple()
                skipped_no_exons += 1

            transcripts.append(
                TxRecord(
                    transcript_id=tid,
                    gene_id=gid,
                    gene_name=gname,
                    seqid=tx.seqid,
                    start=tx.start,
                    end=tx.end,
                    strand=tx.strand,
                    exons=exon_intervals,
                    attrs=attr_dict(tx),
                    gene_type=gtype,
                    transcript_type=ttype,
                    tags=tags,
                )
            )

    log(f"Loaded {len(transcripts)} {label} transcripts; {skipped_no_exons} have no exon children")
    return transcripts


def build_reference_index(ref_txs):
    by_seqid = defaultdict(list)
    for tx in ref_txs:
        by_seqid[tx.seqid].append(tx)
    for seqid in by_seqid:
        by_seqid[seqid].sort(key=lambda tx: tx.start)
    return by_seqid


def candidate_refs(query, ref_index):
    refs = ref_index.get(query.seqid, [])
    for ref in refs:
        if ref.start > query.end:
            break
        if ref.end < query.start:
            continue
        yield ref


def default_reference_hints_output(output):
    if output.endswith(".tsv"):
        return output[:-4] + ".reference_hints.tsv"
    return output + ".reference_hints.tsv"


def default_anchor_loci_output(output):
    if output.endswith(".tsv"):
        return output[:-4] + ".anchor_loci.tsv"
    return output + ".anchor_loci.tsv"


def load_sqanti_classification(path):
    if not path:
        return {}

    log(f"Loading SQANTI classification: {path}")
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return {}

        id_column = None
        for candidate in ("isoform", "isoform_id", "transcript_id", "associated_transcript"):
            if candidate in reader.fieldnames:
                id_column = candidate
                break
        if id_column is None:
            id_column = reader.fieldnames[0]

        out = {}
        for row in reader:
            tid = row.get(id_column, "")
            if tid:
                out[tid] = row

    log(f"Loaded SQANTI rows: {len(out)}")
    return out


def first_present(row, keys, default=""):
    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            return value
    return default


def sqanti_summary(sqanti_rows, transcript_id):
    row = sqanti_rows.get(transcript_id, {})
    return {
        "sqanti_structural_category": first_present(row, ("structural_category", "category")),
        "sqanti_associated_gene": first_present(row, ("associated_gene", "associated_gene_id")),
        "sqanti_associated_transcript": first_present(row, ("associated_transcript", "associated_transcript_id")),
    }


def is_missing(value):
    return value is None or str(value).strip() == ""


def is_novel_transcript_id(value):
    if is_missing(value):
        return True
    return str(value).strip().lower() in {"novel", "na", "nan", "none"}


def reference_hits(tx, ref_index, args, require_opposite_strand=False):
    hits = []
    for ref in candidate_refs(tx, ref_index):
        if ref.strand not in ("+", "-"):
            continue
        if require_opposite_strand and tx.strand == ref.strand:
            continue

        span_bp = span_overlap_bp(tx, ref)
        if span_bp < args.min_span_overlap_bp:
            continue

        exonic_bp = overlap_bp(tx.exons, ref.exons)

        tx_fraction = exonic_bp / tx.exonic_len if tx.exonic_len else 0.0
        ref_fraction = exonic_bp / ref.exonic_len if ref.exonic_len else 0.0
        if tx_fraction < args.min_query_exonic_fraction:
            continue
        if ref_fraction < args.min_reference_exonic_fraction:
            continue

        hits.append(
            {
                "ref": ref,
                "span_overlap_bp": span_bp,
                "exonic_overlap_bp": exonic_bp,
                "tx_exonic_overlap_fraction": tx_fraction,
                "ref_exonic_overlap_fraction": ref_fraction,
            }
        )
    return hits


def best_reference_hit(hits, tx=None, legacy=False):
    """Pick the reference transcript a query transcript most likely belongs to.

    BUG FIX (v1.1.0): the legacy ranking used raw overlapping bp with a
    tie-break on the LONGEST reference span. For a short antisense lncRNA
    embedded in a long host gene, overlap bp ties (the lncRNA is fully covered
    either way) and the tie-break then picks the host gene -- on the opposite
    strand. The pair's two anchors collapse onto the same gene and the antisense
    pair silently disappears at gene level. Exactly the embedded-lncRNA case the
    pipeline exists to find.

    New ranking, in order of priority:
      1. same strand   - a transcript belongs to a gene on its own strand
      2. reciprocal overlap - min(query_fraction, ref_fraction); a genuine match
         covers both partners, containment covers only one
      3. exonic overlap bp
      4. SHORTEST reference span (tighter containment), not the longest

    Pass legacy=True to reproduce pre-1.1.0 anchoring for comparison.
    """
    if not hits:
        return None
    if legacy or tx is None:
        return sorted(
            hits,
            key=lambda hit: (hit["exonic_overlap_bp"], hit["span_overlap_bp"], hit["ref"].span_len),
            reverse=True,
        )[0]

    def rank(hit):
        ref = hit["ref"]
        return (
            1 if ref.strand == tx.strand else 0,
            min(hit["tx_exonic_overlap_fraction"], hit["ref_exonic_overlap_fraction"]),
            hit["exonic_overlap_bp"],
            -ref.span_len,
        )

    return sorted(hits, key=rank, reverse=True)[0]


def anchor_from_sqanti_or_reference(tx, sqanti, hits, legacy_anchoring=False):
    associated_gene = sqanti["sqanti_associated_gene"]
    associated_transcript = sqanti["sqanti_associated_transcript"]
    best_hit = best_reference_hit(hits, tx=tx, legacy=legacy_anchoring)

    if not is_missing(associated_gene) and not is_novel_transcript_id(associated_transcript):
        gene_name = ""
        if best_hit and best_hit["ref"].gene_id == associated_gene:
            gene_name = best_hit["ref"].gene_name
        return {
            "anchor_gene_id": associated_gene,
            "anchor_gene_name": gene_name,
            "anchor_transcript_id": associated_transcript,
            "anchor_source": "sqanti_associated_transcript",
        }

    if not is_missing(associated_gene):
        gene_name = ""
        if best_hit and best_hit["ref"].gene_id == associated_gene:
            gene_name = best_hit["ref"].gene_name
        return {
            "anchor_gene_id": associated_gene,
            "anchor_gene_name": gene_name,
            "anchor_transcript_id": "novel",
            "anchor_source": "sqanti_associated_gene",
        }

    if best_hit:
        ref = best_hit["ref"]
        return {
            "anchor_gene_id": ref.gene_id,
            "anchor_gene_name": ref.gene_name,
            "anchor_transcript_id": ref.transcript_id,
            "anchor_source": "reference_best_hit",
        }

    return {
        "anchor_gene_id": tx.gene_id,
        "anchor_gene_name": tx.gene_name,
        "anchor_transcript_id": "",
        "anchor_source": "query_gtf",
    }


def genes_overlapping_pair_span(tx1, tx2, ref_index):
    start = max(tx1.start, tx2.start)
    end = min(tx1.end, tx2.end)
    if start > end:
        return "", ""

    genes = {}
    for ref in ref_index.get(tx1.seqid, []):
        if ref.start > end:
            break
        if ref.end < start:
            continue
        if overlap_bp(((start, end),), ref.exons) <= 0:
            continue
        if ref.gene_id:
            genes[ref.gene_id] = ref.gene_name

    return ",".join(sorted(genes)), ",".join(genes[gid] for gid in sorted(genes))


def reference_support_class(tx1_ref_hits, tx2_ref_hits):
    if tx1_ref_hits and tx2_ref_hits:
        return "both_transcripts_reference_anchored"
    if tx1_ref_hits or tx2_ref_hits:
        return "one_transcript_reference_anchored"
    return "no_reference_anchor"


def call_reconstructed_pairs(query_txs, ref_index, sqanti_rows, args):
    rows = []
    candidate_pairs = 0
    emitted_pairs = set()
    query_index = build_reference_index(query_txs)

    for tx1 in query_txs:
        if tx1.strand not in ("+", "-"):
            continue

        for tx2 in candidate_refs(tx1, query_index):
            if tx2.transcript_id == tx1.transcript_id:
                continue

            pair_key = tuple(sorted((tx1.transcript_id, tx2.transcript_id)))
            if pair_key in emitted_pairs:
                continue

            if tx2.strand not in ("+", "-"):
                continue
            if tx1.strand == tx2.strand:
                continue

            candidate_pairs += 1
            span_bp = span_overlap_bp(tx1, tx2)
            if span_bp < args.min_span_overlap_bp:
                continue

            all_exon_bp = overlap_bp(tx1.exons, tx2.exons)
            first_exon_bp = overlap_bp(first_exon(tx1), first_exon(tx2))
            last_exon_bp = overlap_bp(last_exon(tx1), last_exon(tx2))
            if first_exon_bp < args.min_first_exon_overlap_bp:
                continue

            emitted_pairs.add(pair_key)
            orientation_class = classify_orientation(tx1, tx2)
            evidence_type = overlap_evidence_type(span_bp, first_exon_bp, last_exon_bp, all_exon_bp)
            tx1_ref_hits = reference_hits(tx1, ref_index, args)
            tx2_ref_hits = reference_hits(tx2, ref_index, args)
            tx1_sqanti = sqanti_summary(sqanti_rows, tx1.transcript_id)
            tx2_sqanti = sqanti_summary(sqanti_rows, tx2.transcript_id)
            tx1_anchor = anchor_from_sqanti_or_reference(
                tx1, tx1_sqanti, tx1_ref_hits, legacy_anchoring=args.legacy_anchoring
            )
            tx2_anchor = anchor_from_sqanti_or_reference(
                tx2, tx2_sqanti, tx2_ref_hits, legacy_anchoring=args.legacy_anchoring
            )
            container_id, contained_id, reciprocal = container_contained(tx1, tx2)
            # DR-01: headline length uses exonic overlap; span only when mature
            # RNAs never touch. The basis is recorded, never silently mixed.
            basis = "exonic" if all_exon_bp > 0 else "span"
            basis_bp = all_exon_bp if all_exon_bp > 0 else span_bp

            rows.append(
                {
                    "pair_id": f"{tx1.transcript_id}|{tx2.transcript_id}",
                    "pair_support": pair_support_class(tx1_anchor, tx2_anchor),
                    "reference_support": reference_support_class(tx1_ref_hits, tx2_ref_hits),
                    "chr": tx1.seqid,
                    "orientation_class": orientation_class,
                    "overlap_evidence_type": evidence_type,
                    "transcript_span_overlap_bp": span_bp,
                    "first_exon_overlap_bp": first_exon_bp,
                    "last_exon_overlap_bp": last_exon_bp,
                    "all_exon_overlap_bp": all_exon_bp,
                    "tx1_transcript_id": tx1.transcript_id,
                    "tx1_gene_id": tx1.gene_id,
                    "tx1_gene_name": tx1.gene_name,
                    "tx1_start": tx1.start,
                    "tx1_end": tx1.end,
                    "tx1_strand": tx1.strand,
                    "tx1_exonic_length": tx1.exonic_len,
                    "tx1_sqanti_structural_category": tx1_sqanti["sqanti_structural_category"],
                    "tx1_sqanti_associated_gene": tx1_sqanti["sqanti_associated_gene"],
                    "tx1_sqanti_associated_transcript": tx1_sqanti["sqanti_associated_transcript"],
                    "tx1_anchor_gene_id": tx1_anchor["anchor_gene_id"],
                    "tx1_anchor_gene_name": tx1_anchor["anchor_gene_name"],
                    "tx1_anchor_transcript_id": tx1_anchor["anchor_transcript_id"],
                    "tx1_anchor_source": tx1_anchor["anchor_source"],
                    "tx2_transcript_id": tx2.transcript_id,
                    "tx2_gene_id": tx2.gene_id,
                    "tx2_gene_name": tx2.gene_name,
                    "tx2_start": tx2.start,
                    "tx2_end": tx2.end,
                    "tx2_strand": tx2.strand,
                    "tx2_exonic_length": tx2.exonic_len,
                    "tx2_sqanti_structural_category": tx2_sqanti["sqanti_structural_category"],
                    "tx2_sqanti_associated_gene": tx2_sqanti["sqanti_associated_gene"],
                    "tx2_sqanti_associated_transcript": tx2_sqanti["sqanti_associated_transcript"],
                    "tx2_anchor_gene_id": tx2_anchor["anchor_gene_id"],
                    "tx2_anchor_gene_name": tx2_anchor["anchor_gene_name"],
                    "tx2_anchor_transcript_id": tx2_anchor["anchor_transcript_id"],
                    "tx2_anchor_source": tx2_anchor["anchor_source"],
                    "tx1_tss": tss(tx1),
                    "tx2_tss": tss(tx2),
                    "tss_distance": abs(tss(tx1) - tss(tx2)),
                    "tx1_tes": tes(tx1),
                    "tx2_tes": tes(tx2),
                    "tes_distance": abs(tes(tx1) - tes(tx2)),
                    # --- L0/curation additions ---------------------------
                    "tx1_gene_type": tx1.gene_type,
                    "tx1_biotype_class": tx1.biotype_class,
                    "tx1_transcript_type": tx1.transcript_type,
                    "tx1_n_exons": tx1.n_exons,
                    "tx1_tags": ";".join(tx1.tags),
                    "tx1_tsl": tx1.tsl,
                    "tx2_gene_type": tx2.gene_type,
                    "tx2_biotype_class": tx2.biotype_class,
                    "tx2_transcript_type": tx2.transcript_type,
                    "tx2_n_exons": tx2.n_exons,
                    "tx2_tags": ";".join(tx2.tags),
                    "tx2_tsl": tx2.tsl,
                    "biotype_pair_class": biotype_pair_class(tx1.biotype_class, tx2.biotype_class),
                    "container_transcript_id": container_id,
                    "contained_transcript_id": contained_id,
                    "embedded_reciprocal": reciprocal,
                    "overlap_length_basis": basis,
                    "overlap_length_bp": basis_bp,
                    "overlap_length_stratum": overlap_length_stratum(basis_bp),
                    "has_exonic_evidence": "1" if all_exon_bp > 0 else "0",
                    "igv_locus": igv_locus(tx1, tx2),
                }
            )

    log(f"Query-vs-query opposite-strand candidate pairs checked: {candidate_pairs}")
    log(f"Reconstructed NAT pairs passing filters: {len(rows)}")
    return rows


def call_reference_hints(query_txs, ref_index, sqanti_rows, reconstructed_pair_ids, args):
    rows = []
    candidate_pairs = 0

    for query in query_txs:
        if query.strand not in ("+", "-"):
            continue

        sqanti = sqanti_summary(sqanti_rows, query.transcript_id)

        for hit in reference_hits(query, ref_index, args, require_opposite_strand=True):
            ref = hit["ref"]

            candidate_pairs += 1
            span_bp = hit["span_overlap_bp"]
            exonic_bp = hit["exonic_overlap_bp"]
            first_exon_bp = overlap_bp(first_exon(query), first_exon(ref))
            last_exon_bp = overlap_bp(last_exon(query), last_exon(ref))
            query_exonic_len = query.exonic_len
            ref_exonic_len = ref.exonic_len
            query_exonic_fraction = hit["tx_exonic_overlap_fraction"]
            ref_exonic_fraction = hit["ref_exonic_overlap_fraction"]

            orientation_class = classify_orientation(query, ref)
            evidence_type = overlap_evidence_type(span_bp, first_exon_bp, last_exon_bp, exonic_bp)

            rows.append(
                {
                    "hint_support": (
                        "query_has_reconstructed_nat_partner"
                        if query.transcript_id in reconstructed_pair_ids
                        else "reference_only_antisense_hint"
                    ),
                    "query_transcript_id": query.transcript_id,
                    "query_gene_id": query.gene_id,
                    "query_gene_name": query.gene_name,
                    "query_chr": query.seqid,
                    "query_start": query.start,
                    "query_end": query.end,
                    "query_strand": query.strand,
                    "query_exonic_length": query_exonic_len,
                    "query_sqanti_structural_category": sqanti["sqanti_structural_category"],
                    "query_sqanti_associated_gene": sqanti["sqanti_associated_gene"],
                    "query_sqanti_associated_transcript": sqanti["sqanti_associated_transcript"],
                    "ref_transcript_id": ref.transcript_id,
                    "ref_gene_id": ref.gene_id,
                    "ref_gene_name": ref.gene_name,
                    "ref_chr": ref.seqid,
                    "ref_start": ref.start,
                    "ref_end": ref.end,
                    "ref_strand": ref.strand,
                    "ref_exonic_length": ref_exonic_len,
                    "orientation_class": orientation_class,
                    "overlap_evidence_type": evidence_type,
                    "transcript_span_overlap_bp": span_bp,
                    "first_exon_overlap_bp": first_exon_bp,
                    "last_exon_overlap_bp": last_exon_bp,
                    "all_exon_overlap_bp": exonic_bp,
                    "query_exonic_overlap_fraction": f"{query_exonic_fraction:.6f}",
                    "ref_exonic_overlap_fraction": f"{ref_exonic_fraction:.6f}",
                    "query_tss": tss(query),
                    "ref_tss": tss(ref),
                    "tss_distance": abs(tss(query) - tss(ref)),
                    "query_tes": tes(query),
                    "ref_tes": tes(ref),
                    "tes_distance": abs(tes(query) - tes(ref)),
                }
            )

    log(f"Query-vs-reference opposite-strand hints checked: {candidate_pairs}")
    log(f"Reference antisense hints passing filters: {len(rows)}")
    return rows


PAIR_FIELDNAMES = [
    "pair_id",
    "pair_support",
    "reference_support",
    "chr",
    "orientation_class",
    "overlap_evidence_type",
    "transcript_span_overlap_bp",
    "first_exon_overlap_bp",
    "last_exon_overlap_bp",
    "all_exon_overlap_bp",
    "tx1_transcript_id",
    "tx1_gene_id",
    "tx1_gene_name",
    "tx1_start",
    "tx1_end",
    "tx1_strand",
    "tx1_exonic_length",
    "tx1_sqanti_structural_category",
    "tx1_sqanti_associated_gene",
    "tx1_sqanti_associated_transcript",
    "tx1_anchor_gene_id",
    "tx1_anchor_gene_name",
    "tx1_anchor_transcript_id",
    "tx1_anchor_source",
    "tx2_transcript_id",
    "tx2_gene_id",
    "tx2_gene_name",
    "tx2_start",
    "tx2_end",
    "tx2_strand",
    "tx2_exonic_length",
    "tx2_sqanti_structural_category",
    "tx2_sqanti_associated_gene",
    "tx2_sqanti_associated_transcript",
    "tx2_anchor_gene_id",
    "tx2_anchor_gene_name",
    "tx2_anchor_transcript_id",
    "tx2_anchor_source",
    "tx1_tss",
    "tx2_tss",
    "tss_distance",
    "tx1_tes",
    "tx2_tes",
    "tes_distance",
    # --- L0/curation additions, appended so column order stays stable -----
    "tx1_gene_type",
    "tx1_biotype_class",
    "tx1_transcript_type",
    "tx1_n_exons",
    "tx1_tags",
    "tx1_tsl",
    "tx2_gene_type",
    "tx2_biotype_class",
    "tx2_transcript_type",
    "tx2_n_exons",
    "tx2_tags",
    "tx2_tsl",
    "biotype_pair_class",
    "container_transcript_id",
    "contained_transcript_id",
    "embedded_reciprocal",
    "overlap_length_basis",
    "overlap_length_bp",
    "overlap_length_stratum",
    "has_exonic_evidence",
    "igv_locus",
    "anchor_locus_id",
    "n_anchor_genes_in_locus",
    "n_pairs_in_anchor_locus",
]


HINT_FIELDNAMES = [
        "hint_support",
        "query_transcript_id",
        "query_gene_id",
        "query_gene_name",
        "query_chr",
        "query_start",
        "query_end",
        "query_strand",
        "query_exonic_length",
        "query_sqanti_structural_category",
        "query_sqanti_associated_gene",
        "query_sqanti_associated_transcript",
        "ref_transcript_id",
        "ref_gene_id",
        "ref_gene_name",
        "ref_chr",
        "ref_start",
        "ref_end",
        "ref_strand",
        "ref_exonic_length",
        "orientation_class",
        "overlap_evidence_type",
        "transcript_span_overlap_bp",
        "first_exon_overlap_bp",
        "last_exon_overlap_bp",
        "all_exon_overlap_bp",
        "query_exonic_overlap_fraction",
        "ref_exonic_overlap_fraction",
        "query_tss",
        "ref_tss",
        "tss_distance",
        "query_tes",
        "ref_tes",
        "tes_distance",
]


def write_manifest(args, pair_rows, hint_rows, hints_output, runtime_sec):
    """Provenance sidecar: without it a truth set cannot be reproduced."""
    manifest = {
        "tool": "nat_caller.py",
        "version": __version__,
        "argv": sys.argv,
        "params": {
            "query": os.path.abspath(args.query),
            "reference": os.path.abspath(args.reference),
            "sqanti_classification": args.sqanti_classification,
            "min_span_overlap_bp": args.min_span_overlap_bp,
            "min_first_exon_overlap_bp": args.min_first_exon_overlap_bp,
            "min_query_exonic_fraction": args.min_query_exonic_fraction,
            "min_reference_exonic_fraction": args.min_reference_exonic_fraction,
            "legacy_anchoring": args.legacy_anchoring,
        },
        "outputs": {
            "pairs": os.path.abspath(args.output),
            "reference_hints": os.path.abspath(hints_output),
        },
        "counts": {"pairs": len(pair_rows), "reference_hints": len(hint_rows)},
        "runtime_sec": round(runtime_sec, 2),
    }
    path = args.output + ".manifest.json"
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    return path


def assign_anchor_loci(pair_rows):
    """Group reconstructed pairs into loci by ANCHOR gene, not raw gene_id.

    Raw query gene_id is often a throwaway FLAIR-assigned identifier for de
    novo transcripts and is not a stable grouping key. anchor_gene_id is the
    same field nat_curate.py's assign_loci() groups by conceptually, just
    computed from SQANTI/reference anchoring instead of a static annotation.
    Same union-find idea, kept intentionally small: no new dependency, no new
    output format beyond one companion TSV.
    """
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx

    for row in pair_rows:
        union(row["tx1_anchor_gene_id"], row["tx2_anchor_gene_id"])

    members = defaultdict(set)
    for gene_id in list(parent):
        members[find(gene_id)].add(gene_id)

    locus_of_anchor = {}
    locus_info = {}
    for idx, (root, gene_ids) in enumerate(sorted(members.items(), key=lambda kv: sorted(kv[1])), 1):
        locus_id = f"ALOC{idx:04d}"
        locus_info[locus_id] = sorted(gene_ids)
        for gene_id in gene_ids:
            locus_of_anchor[gene_id] = locus_id

    pairs_per_locus = Counter()
    for row in pair_rows:
        locus_id = locus_of_anchor[row["tx1_anchor_gene_id"]]
        row["anchor_locus_id"] = locus_id
        pairs_per_locus[locus_id] += 1

    for row in pair_rows:
        row["n_pairs_in_anchor_locus"] = pairs_per_locus[row["anchor_locus_id"]]
        row["n_anchor_genes_in_locus"] = len(locus_info[row["anchor_locus_id"]])

    locus_rows = [
        {
            "locus_id": locus_id,
            "n_anchor_genes": len(gene_ids),
            "n_pairs_in_locus": pairs_per_locus.get(locus_id, 0),
            "complex_locus": "1" if len(gene_ids) > 2 else "0",
            "anchor_gene_ids": ",".join(gene_ids),
        }
        for locus_id, gene_ids in sorted(locus_info.items())
    ]
    n_complex = sum(1 for r in locus_rows if r["complex_locus"] == "1")
    log(f"Anchor loci: {len(locus_rows)} connected components, {n_complex} with >2 anchor genes")
    return locus_rows


def write_tsv(rows, output, fieldnames):
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    start_time = time.time()
    db_paths = []

    try:
        query_db_path, query_db = create_feature_db(
            args.query,
            args.output,
            "query",
            tmpdir=args.tmpdir,
            keep_db=args.keep_db,
        )
        ref_db_path, ref_db = create_feature_db(
            args.reference,
            args.output,
            "reference",
            tmpdir=args.tmpdir,
            keep_db=args.keep_db,
        )
        db_paths.extend([query_db_path, ref_db_path])

        query_txs = load_transcripts(query_db, "query")
        ref_txs = load_transcripts(ref_db, "reference")
        ref_index = build_reference_index(ref_txs)
        sqanti_rows = load_sqanti_classification(args.sqanti_classification)

        pair_rows = call_reconstructed_pairs(query_txs, ref_index, sqanti_rows, args)
        anchor_locus_rows = assign_anchor_loci(pair_rows)
        write_tsv(pair_rows, args.output, PAIR_FIELDNAMES)
        anchor_loci_output = default_anchor_loci_output(args.output)
        write_tsv(
            anchor_locus_rows, anchor_loci_output,
            ["locus_id", "n_anchor_genes", "n_pairs_in_locus", "complex_locus", "anchor_gene_ids"],
        )

        reconstructed_pair_ids = set()
        for row in pair_rows:
            reconstructed_pair_ids.add(row["tx1_transcript_id"])
            reconstructed_pair_ids.add(row["tx2_transcript_id"])

        hints_output = args.reference_hints_output or default_reference_hints_output(args.output)
        hint_rows = call_reference_hints(query_txs, ref_index, sqanti_rows, reconstructed_pair_ids, args)
        write_tsv(hint_rows, hints_output, HINT_FIELDNAMES)

        manifest_path = write_manifest(
            args, pair_rows, hint_rows, hints_output, time.time() - start_time
        )

        log(f"Wrote: {args.output}")
        log(f"Wrote: {hints_output}")
        log(f"Wrote: {anchor_loci_output}")
        log(f"Wrote: {manifest_path}")
        log(f"Total runtime: {time.time() - start_time:.2f} sec")

    finally:
        if not args.keep_db:
            for db_path in db_paths:
                try:
                    os.remove(db_path)
                except FileNotFoundError:
                    pass
                except sqlite3.Error:
                    pass


if __name__ == "__main__":
    main()