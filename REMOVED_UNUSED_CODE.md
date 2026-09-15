# Disconnected / unused / stub code

Audit performed 2026-09-15 by walking every module import against its actual
invocation site, every `params.*` declaration against its actual read sites, and
every named profile's samplesheet against the current input schema — not from
memory, each item below was independently re-confirmed against the live code
before touching anything.

**What "removed" means here**: nothing was deleted from disk. Everything below
was disconnected from the live workflow graph (its `include {}` line, its
`nextflow.config` block, and/or its `params.*` declaration removed) and, where
it's a whole file/directory, added to `.gitignore` so it doesn't get published —
but it's all still sitting in the working tree in case it's worth finishing or
reviving later. This file is the map of what's there and why it's inert.

---

## 1. Orphaned modules (imported, never invoked)

### `SEQKIT_CONCAT` (modules/local/seqkit/concat/)

Was imported in `subworkflows/local/preprocess_long_reads/main.nf` but never
called — `FASTQ_CONCAT` (a plain `cat`) is what actually reassembles chunked
fastqs there. Had its own `withName: SEQKIT_CONCAT {}` resource/publishDir
block in `nextflow.config`, also unused. Done:
- removed the `include { SEQKIT_CONCAT }` line from `preprocess_long_reads/main.nf`
- removed its `nextflow.config` block (was: `publishDir` → `fastq/trimmed`,
  `symlink`, `ext.prefix = "${meta.id}_concat"`)
- `modules/local/seqkit/concat/` itself (`main.nf`, `environment.yml`,
  `meta.yml`, `tests/`) — left on disk, added to `.gitignore`

(Note: `modules/nf-core/seqkit/concat/` — the separate, unmodified nf-core copy
this local one was based on — was left completely alone, including in
`modules.json`. It's a different directory, also not imported by anything, but
touching an nf-core-vendored module felt like a separate decision from
disconnecting this repo's own local copy.)

### `CREATE_MANIFESTS` (modules/local/flair/create_manifests/)

Was imported in `subworkflows/local/flair_combine_quantify/main.nf` but never
called — that subworkflow builds its FLAIR manifests inline in each process's
own bash/`exec:` block instead (see `modules/local/flair/combine/main.nf` and
`modules/local/flair/quantify/main.nf`). Looks like an earlier, alternative
implementation of the same manifest-building step (pure-Groovy `exec:` block
using `toRealPath()`) that was superseded but never deleted. Done:
- removed the `include { CREATE_MANIFESTS }` line from
  `flair_combine_quantify/main.nf`
- removed its `nextflow.config` block (was: `publishDir` →
  `flair/MANIFESTS/${meta}`, `link`)
- `modules/local/flair/create_manifests/` itself — left on disk, added to
  `.gitignore`

### `SORTMERNA` (modules/nf-core/sortmerna/)

Was imported in `subworkflows/local/preprocess_short_reads/main.nf`, but its
call site was an empty, commented-out block:
```groovy
/*
if (!params.skip_sortmerna) {
    
}
*/
```
i.e. the `--skip_sortmerna` param existed and was checked, but the branch it
guarded never did anything — SortMeRNA has never actually run in this
pipeline. Done:
- removed the `include { SORTMERNA }` line and the empty commented-out block
  above, from `preprocess_short_reads/main.nf`
- removed `params.skip_sortmerna` from `nextflow.config`
- `modules/nf-core/sortmerna/` itself and its `modules.json` entry — directory
  left on disk (added to `.gitignore`); the `modules.json` entry *was* removed
  outright (unlike the directory, `modules.json` is published, and leaving it
  claim a module is installed when the directory won't be there for anyone who
  clones the repo would be actively misleading, not just inert)

### `MULTIQC` (modules/nf-core/multiqc/)

Installed (`nf-core modules install`) but never imported by any workflow or
subworkflow — confirmed via a repo-wide sweep, not just checked near its own
directory. `preprocess_short_reads/main.nf` did build a `ch_multiqc_files`
channel (mixing `FASTQC_RAW`/`FASTP`/`TRIMMOMATIC`/`FASTQC_TRIMMED` outputs) and
emit it as `multiqc_files` — clearly scaffolding for a MultiQC step that was
planned but never wired up: that emit was never read by `workflows/natseq.nf`
either. The `multiqc_report.html`/`multiqc_data/` that used to sit at the repo
root were from a one-off manual `multiqc` CLI invocation, not this pipeline.
Done:
- removed `ch_multiqc_files` construction and the `multiqc_files` emit from
  `preprocess_short_reads/main.nf` (kept `FASTQC_RAW`/`FASTP`/`TRIMMOMATIC`/
  `FASTQC_TRIMMED` themselves — those do run and are useful on their own, only
  the dead multiqc-file-collection plumbing around them was removed)
- `modules/nf-core/multiqc/` directory left on disk, added to `.gitignore`;
  its `modules.json` entry removed outright (same reasoning as SORTMERNA above)

---

## 2. Params declared but never read anywhere

Confirmed by grepping the entire repo (`.nf` + `.config`) for `params.<name>` —
zero hits outside the declaration itself. All four removed from `nextflow.config`:

| Param | Was | Why it's gone |
|---|---|---|
| `skip_flair_comb_quant` | `false` | Declared for skipping FLAIR combine/quantify; `workflows/natseq.nf` calls `FLAIR_COMBINE_QUANTIFY(...)` unconditionally — this flag was never actually checked anywhere |
| `skip_sqanti` | `false` | Same story for `SQANTI(...)` — always runs, flag never checked |
| `skip_sortmerna` | `false` | See SORTMERNA above — guarded an empty block |
| `f_quantify_args` | `''` | Intended as free-form extra args for `flair quantify`, matching the pattern used by `f_combine_args`/`filtlong_args`/etc. But `modules/local/flair/quantify/main.nf` reads `task.ext.args`, and `nextflow.config`'s `withName: FLAIR_QUANTIFY {}` block only sets `publishDir` — never `ext.args`. So this param was declared but had no path into the actual `flair quantify` command. **This one would be a one-line fix rather than a real dead end** if you want it working instead of removed — add to the `withName: FLAIR_QUANTIFY` block: `ext.args = { params.f_quantify_args }`. |

---

## 3. Broken/non-portable Nextflow profiles

All six named profiles that used to live in `nextflow.config`'s `profiles {}`
block are gone (the block itself is gone too, replaced with a one-line comment
pointing here). Two different reasons:

**`test_diab`, `test_gold`, `test_sgnex`** — genuinely broken, not just
stale paths. They pointed at `TEST/SAMPLESHEETS/*.csv`, whose header schemas
predate the current input schema entirely:
```
samplesheet_test_DIABETES_3_S.csv:  sample,type,fastq
samplesheet_test_GOLD.csv:          sample,tissue_type,read_type,fastq_1,fastq_2
samplesheet_test_SGNEX.csv:         sample,tissue_type,read_type,fastq_1,fastq_2
```
`workflows/natseq.nf` reads `row.sample_id`, `row.group_id`, `row.library_type`,
etc. — none of which exist in these files. Running any of these three profiles
as-is fails immediately (`library_type` resolves to `null`, which then throws in
the `MINIMAP2_ALIGN` preset selector).

**`diab`, `deme`, `sgnex`** — pointed at root-level samplesheets
(`samplesheet.csv`, `samplesheet_DTM.csv`, `samplesheet_SGNEX.csv`) whose every
row hardcodes `/home/artale/Work/NATs/...` — a different machine, not this one,
not present anywhere in this repo or its `.gitignore`d data. `samplesheet_DTM.csv`
additionally uses lowercase `library_type` values (`cdna`) that don't
case-sensitively match the pipeline's own `dRNA`/`cDNA`/`direct-cDNA` checks, so
even a path fix alone wouldn't make it runnable.

Removed the entire `profiles {}` block:
```groovy
profiles {
    test_diab {
        params.input = "${projectDir}/TEST/SAMPLESHEETS/samplesheet_test_DIABETES_3_S.csv"
        process.maxForks = 6
    }

    test_gold {
        params.input = "${projectDir}/TEST/SAMPLESHEETS/samplesheet_test_GOLD.csv"
    }

    test_sgnex {
        params.input = "${projectDir}/TEST/SAMPLESHEETS/samplesheet_test_SGNEX.csv"
    }


    diab {
        params.input = "${projectDir}/samplesheet.csv"
    }

    deme {
        params.input = "${projectDir}/samplesheet_DTM.csv"
    }

    sgnex {
        params.input = "${projectDir}/samplesheet_SGNEX.csv"
    }

}
```
The underlying files (`samplesheet*.csv`, `TEST/`) were left on disk
untouched (and were already excluded from git before this pass) — only the
`nextflow.config` profile blocks that pointed at them are gone. Running the
pipeline now means always passing `--input` explicitly (which is what every
actual run command used in practice anyway — see `README.md`).

---

## 4. Dead / broken scripts in `bin/`

### `bin/script_extract_overlap.py`

Syntactically invalid Python — `filename = ` and `dbfname = ` (lines 9-10) are
incomplete assignments, so the file doesn't even pass parsing. Never called by
any process (only `..._mod.py`, and later `nat_caller.py`, ever were). Looks
like the very first draft that `script_extract_overlap_mod.py` was written
from. Left on disk, added to `.gitignore`.

### `bin/script_extract_overlap_mod.py`

Was called from `FIND_OVERLAP` earlier in this project's history (produced
`${meta.id}.overlaps.csv`, a cruder overlap list restricted to *adjacent* genes
by sort order — see the module's difference from `nat_caller.py`'s full
candidate scan). `modules/local/find_overlap/main.nf` now calls only
`nat_caller.py`; this script is no longer invoked from anywhere. Left on disk,
added to `.gitignore`.

### `nat_caller_summary.R`

Never wired into Nextflow — a standalone, manually-run R script. Its hardcoded
`input_files`/`outdir` paths (`~/Desktop/MAGA/DEMENTIA/...`) point at an
unrelated project/directory on the original author's machine, not this
repository's data. Already excluded from git via `.gitignore` before this pass
(no change needed here). If you want a "run this after `nat_caller.py`
completes" summary/plotting script for this repo specifically, it would need to
be rewritten with CLI args instead of hardcoded paths — this version isn't a
useful starting point beyond its plotting logic.

---

## 5. Left in place, flagged only (nothing changed)

- **`bin/intron-prospector-merge`** and its wrapper **`bin/intronProspectorMerge`**
  — confirmed unused (not called by any `.nf` process, unlike
  `bin/intron-prospector` / `bin/intronProspector` which the `INTRONPROSPECTOR`
  module does call). Already `.gitignore`d, not touched further this pass.
- **`null.overlaps.csv`** (repo root) — a bug artifact from a run where
  `meta.id` resolved to the literal string `"null"`. Not code, already
  `.gitignore`d; left on disk as-is since this pass was scoped to code, not
  leftover run output.
