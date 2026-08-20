# annotate_snps

One cohort VCF -> CADD, SpliceAI, SQUIRLS and ANNOVAR annotations, merged to a single
multianno TSV.

Input is an already-merged VCF, not a samplesheet -- this pipeline annotates what it is
given. To annotate more than one variant set, run it again with a separate `--vcf`/`--outdir`
per set.

Chromosome intervals are derived at runtime from what is actually in the input VCF
(`list_chroms`), not from a static `chr1..chrY` list: on targeted-panel data most
chromosomes have zero variants, and CADD/SpliceAI/SQUIRLS all error out on an empty VCF.

## Run

**Needs Nextflow 26.04.0 or newer.** Typed syntax (`nextflow.enable.types`) is a preview
feature that does not exist in 25.10.x, and the manifest's `nextflowVersion` will not tell you
so: Nextflow compiles the modules before it checks that field, so an older release dies with
`Unrecognized feature flag 'nextflow.enable.types'` followed by every include reported as
"Module could not be parsed". Measured against 25.10.0 (fails) and 26.04.6 (passes).

Export the v2 parser before running:

```bash
export NXF_SYNTAX_PARSER=v2
```

```bash
nextflow run nf/annotate_snps/main.nf \
    -profile <mygroup>,compute2 \
    --vcf cohort.vcf.gz \
    --tbi cohort.vcf.gz.tbi \
    --cohort WURARE \
    --data_type wgs \
    --reference.fa /path/to/the_reference_the_vcf_was_called_against.fa \
    --cadd_data_dir /storage1/fs1/WURare/Active/annotations/CADD/v1.6/ \
    --squirls_config /storage1/fs1/WURare/Active/annotations/squirls/c2.squirls-config.yml \
    --squirls_jannovar_model /storage1/fs1/WURare/Active/annotations/squirls/jannovar_v0.35/hg38_refseq.ser \
    --spliceai_precomputed_scores /storage1/fs1/WURare/Active/annotations/spliceai/splice.vcf.gz \
    --spliceai_precomputed_tbi /storage1/fs1/WURare/Active/annotations/spliceai/splice.vcf.gz.tbi \
    --annovar_dir /storage1/fs1/WURare/Active/annotations/annovar/annovar/ \
    --annovar_splice_scores_script /path/to/annovar_parse_out_spliceai_squirls_scores.py \
    --omim_xref /storage1/fs1/WURare/Active/annotations/omim/20260702/gene_xref.annovar.txt \
    --outdir results
```

Selecting no profile gives a local, container-less run. `-profile compute2` enables apptainer
and SLURM submission; the group profile supplies the SLURM account -- use `wurare` for this repo's runs
(`cooper` charges a different account). Cluster setup is self-contained -- this repo ships the
profile at [`../shared/conf/compute2.config`](../shared/conf/compute2.config).

Every required param with no default is checked by name at startup, so a missing one fails
immediately rather than surfacing later as a `checkIfExists` on a path you never set.

## Params

Required, no default -- the run fails at startup without them. Five are required only when the
stage that reads them runs; see the "Required when" column and the re-annotation section:

| Param | What it is | Required when |
|---|---|---|
| `--vcf` / `--tbi` | Cohort VCF to annotate and its tabix index | always |
| `--cohort` / `--data_type` | Output filename components, e.g. `WURARE` / `wgs` | always |
| `--annovar_dir` | ANNOVAR install, containing `humandb/` and the perl scripts | always |
| `--annovar_splice_scores_script` | `annovar_parse_out_spliceai_squirls_scores.py` | always |
| `--omim_xref` | OMIM gene xref table for `table_annovar.pl -xref` | always |
| `--cadd_data_dir` | CADD v1.6 GRCh38 data directory, must contain `CADD.sh` | unless `--precomputed_cadd` **and** `--skip_spliceai_squirls` |
| `--squirls_config` | SQUIRLS config YAML | unless `--skip_spliceai_squirls` |
| `--squirls_jannovar_model` | SQUIRLS Jannovar transcript model (`.ser`) | unless `--skip_spliceai_squirls` |
| `--spliceai_precomputed_scores` | Precomputed SpliceAI scores VCF | unless `--skip_spliceai_squirls` |
| `--spliceai_precomputed_tbi` | Its tabix index | unless `--skip_spliceai_squirls` |

Defaulted:

| Param | Default | Notes |
|---|---|---|
| `--outdir` | `results` | |
| `--clinvar_date` | `clinvar_20260627` | ANNOVAR protocol name **and** the source of the run date in every output filename -- `annovar.nf` splits on `_` and takes the second token, giving `<cohort>_<data_type>_20260627`. Must name a database actually present in `${annovar_dir}/humandb/` |
| `--annovar_protocols` / `--annovar_operations` | 16 protocols; see `conf/params.config` | Paired **position for position**, so the two must stay the same length -- `main.nf` refuses to start otherwise, because a mismatch silently applies the wrong operation to every database after the gap. Overriding `protocols` means overriding `operations` too. `CLINVAR` in the list is a placeholder substituted with `--clinvar_date`, and every other entry must exist in `humandb/` (the startup preflight checks). Detail under [Preparing external inputs](#preparing-external-inputs) below |
| `--reference.fa` / `.fai` / `.dict` | iGenomes GATK.GRCh38 (`s3://`) | See the warning below |
| `--somalier_sites` | `null` | Opt-in sample QC. GRCh38 sites VCF: [sites.hg38.vcf.gz](https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz), listed with the [releases](https://github.com/brentp/somalier/releases). This alone is enough to run the QC -- `--ped` is optional |
| `--ped` | `null` | Standard 6-column PED, and **optional**. With one, measured relatedness and sex are checked *against* the pedigree. Without one, somalier infers the family structure (`somalier relate --infer`) and that inference is checked instead. Either way `04_qc` gains `ploidy.tsv` and `cohort.inferred.ped`. Unused unless `--somalier_sites` is set |
| `--somalier_labels` / `--somalier_1kg_dir` | `null` | Optional ancestry, **both or neither**, and both need `--somalier_sites`. [ancestry-labels-1kg.tsv](https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv) and the unpacked [1kg.somalier.tar.gz](https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz) — the param names the **directory**, not the tarball |
| `--alignments` | `null` | Optional CSV `sample,alignment,alignment_index` (CRAM or BAM, same sheet as `annotate_svs`). When set, somalier extracts depth from the alignments instead of `FORMAT/AD` in `--vcf` — use it for exomes, where the filtered VCF under-calls sex. Needs `--somalier_sites` |
| `--alignment_reference` / `--alignment_reference_index` | `null` | Required with `--alignments`: the exact FASTA (+ .fai) the alignments were written against — deliberately separate from `--reference`, CRAM decodes against it |

`--reference.fa` is not authoritative for any particular dataset. `normalize_vcf` runs
`bcftools norm --check-ref s`, which **rewrites REF against whatever this points at** rather
than failing on a mismatch. Override it to the reference your input VCF was actually called
against. The s3 default also means re-fetching per task; download once and pass a local path.

## Sample QC (optional)

`--somalier_sites` turns on a stage that has nothing to do with annotation: it genotypes
~17k common sites out of `--vcf` and reports **observed relatedness** for every pair of
samples, and optionally ancestry. Nothing downstream consumes it — this pipeline has no depth
filter to make ploidy-aware, unlike `annotate_svs` — so it is published and nothing more. It is
the sample-swap check this pipeline otherwise has none of, and a swapped sample is invisible in
an annotated VCF.

**`--ped` is optional**, and what it changes is what the measurements are checked *against*.
With one, somalier is told what the pedigree claims, and the stated structure is checked —
never invented. Without one, somalier infers the family structure itself (`somalier relate
--infer`); sex is still derived from depth by the same check either way. Only where depth alone
cannot decide does the published sex fall back to what the pedigree already carried — the
operator's own value with `--ped`, or somalier's own (non-depth) inference without one. Both
routes publish:

- `04_qc/ploidy.tsv` — per-sample X/Y copies against the pedigree's sex column, stated or
  inferred. Where the two disagree, **the data wins**: the row is marked `DISAGREES` and
  warned about on stderr.
- `04_qc/cohort.inferred.ped` — the pedigree with column 5 replaced by the measured sex.
  With `--ped` the family structure is the operator's, untouched; without one it is
  somalier's inference. Either way, this is the file to hand `annotate_svs` as its `--ped`
  when a cohort arrives with no recorded pedigree.

By default, depth comes from `FORMAT/AD` in the VCF rather than from reads — and on a
**filtered callset that is not enough**. An exome VCF carries too few informative X/Y sites,
every X depth ratio drifts toward 1.0 and every Y toward 0, and the whole cohort is reported
female. For exomes (or any filtered callset), pass the alignments instead:

```
--alignments alignments.csv \
--alignment_reference /path/to/reference.fa \
--alignment_reference_index /path/to/reference.fa.fai
```

`alignments.csv` is the same sheet `annotate_svs` uses — `sample,alignment,alignment_index`,
CRAM or BAM:

```
sample,alignment,alignment_index
FAM01-01,/path/FAM01-01.cram,/path/FAM01-01.cram.crai
```

With a sheet, somalier reads real depth off the alignments (one task per sample) and only
the extract step changes — relatedness, ancestry and the report are fed the same files.
`--alignment_reference` must be the **exact FASTA the alignments were written against**
(CRAM decodes against it; a mismatch gives wrong depth, not an error), which is why it is a
separate knob from `--reference`. The sheet's sample names must match both the alignments'
`@RG SM` (enforced per task) and the VCF header (so the per-sample tables line up). All
three params travel together and `--alignments` also requires `--somalier_sites`; the run
is rejected at startup otherwise.

Without a sheet, a VCF whose caller did not write `AD` genotypes to zero depth, which fails
the run with that named as the likely cause rather than reporting a cohort of undetermined
samples. Relatedness needs only the genotypes and is unaffected either way.

The same processes and the same `check_somalier.awk` serve `annotate_svs` — they live in
`nf/shared/`, so the two pipelines cannot drift.

## Preparing external inputs

Three inputs have no public source and no default. They are host-local files on compute2,
named here rather than hidden, on the principle of stating what is
undocumented:

| Input | Status |
|---|---|
| `--annovar_dir` | A licensed ANNOVAR install plus a populated `humandb/`. Not redistributable; register at the ANNOVAR site for the perl scripts, then download the databases named in the `table_annovar.pl` protocol list |
| `--annovar_splice_scores_script` | `annovar_parse_out_spliceai_squirls_scores.py`, host-local under the operator's scripts directory. Not vendored in any repo and has no upstream. `file()` accepts an `https://` URL pinned to a commit SHA as well as a local path, so this can be pointed at a raw file in git once it lives somewhere |
| `--omim_xref` | OMIM gene xref table. OMIM requires a license; the file is generated per OMIM release (the compute2 copy is dated `20260702`) |

**The ANNOVAR databases are preflighted at startup**, so a missing one fails in seconds rather
than after hours of CADD and SpliceAI. The run stops before submitting a single task and names
every missing database at once, so they can be fixed in one pass.

The protocol list lives in `conf/params.config` as `--annovar_protocols`, paired position for
position with `--annovar_operations`. It is one list, read both by the preflight and by
`table_annovar.pl`, so the check cannot drift from the command. `CLINVAR` in that list is a
placeholder substituted with `--clinvar_date`, which is why overriding the ClinVar release also
moves which database is required.

Three things the check knows that a plain `ls hg38_*` does not:

- **`CADDv1.6` is not expected beforehand** -- `build_cadd_humandb` creates it at runtime and it
  is symlinked in, so the preflight skips it.
- **`1000g2015aug_all` does not follow the naming convention.** ANNOVAR maps it onto
  `hg38_ALL.sites.2015_08.txt`; a `hg38_1000g2015aug_all.txt` check would fail on a correctly
  populated install.
- **`refGene` needs `hg38_refGeneMrna.fa`** as well as its `.txt`, because its operation is `gx`.

`.idx` files are deliberately **not** required per protocol -- ANNOVAR falls back to a linear
scan without one, so demanding them would be stricter than the tool itself. At least one `.idx`
must exist somewhere, because `stage_humandb` globs `*.txt.idx` and `ln -s` fails on an
unmatched glob. The four perl scripts are checked too.

To see the same picture by hand:

```bash
ls ${annovar_dir}/humandb/hg38_*
```

## Output

```
<outdir>/
  01_cadd/
    caddv1.6.out.tsv.gz            merged genome-wide CADD table
    caddv1.6.out.tsv.gz.tbi
  02_annovar/
    <cohort>_<data_type>_<date>.hg38_multianno.tsv
  03_annotated_vcf/
    <cohort>_<data_type>.spliceai_squirls.vcf.gz    SpliceAI + SQUIRLS INFO only
    <cohort>_<data_type>.spliceai_squirls.vcf.gz.tbi
  04_qc/                             only with --somalier_sites
    ploidy.tsv                       X/Y copies, pedigree sex (stated or inferred), agreement
    cohort.inferred.ped              column 5 replaced by the measured sex; family structure
                                     is the operator's PED, or somalier's inference without one
    <cohort>.samples.tsv             somalier's own tables
    <cohort>.pairs.tsv               pairwise relatedness; the expected column needs --ped
    <cohort>.html
    somalier-ancestry.tsv            only with --somalier_labels/--somalier_1kg_dir;
                                     this run's samples AND the 1kg reference set
    <cohort>.somalier-ancestry.tsv   the same table, this run's samples only
    somalier-ancestry.html           the PCA plot
  pipeline_info/
    software_versions.yml          one entry per process, deduplicated
```

`03_annotated_vcf` is numbered after `02_annovar` but produced before it -- renumbering the
existing directories would break paths already documented here.

### Re-annotating without a full re-run

`03_annotated_vcf` exists so ANNOVAR can be re-run against updated databases without
repeating CADD, SpliceAI and SQUIRLS, which are the expensive part of a run. Together with
`01_cadd/caddv1.6.out.tsv.gz` it is everything the ANNOVAR stage consumes.

**Read the filename literally: this VCF carries SpliceAI and SQUIRLS INFO fields and nothing
else.** CADD scores are not in it -- CADD is emitted as its own TSV and joined by ANNOVAR at
the `humandb/` level, and unlike upstream this pipeline has no vcfanno step to fold
CADD into INFO. No ANNOVAR field is in it either; `table_annovar.pl` writes TSV, never VCF.

To re-run ANNOVAR alone against updated databases, point a new run at both artifacts:

```bash
nextflow run nf/annotate_snps/main.nf \
    -profile <mygroup>,compute2 \
    --vcf   results/03_annotated_vcf/WURARE_wgs.spliceai_squirls.vcf.gz \
    --tbi   results/03_annotated_vcf/WURARE_wgs.spliceai_squirls.vcf.gz.tbi \
    --skip_spliceai_squirls true \
    --precomputed_cadd     results/01_cadd/caddv1.6.out.tsv.gz \
    --precomputed_cadd_tbi results/01_cadd/caddv1.6.out.tsv.gz.tbi \
    --clinvar_date clinvar_20270101 \
    --cohort WURARE --data_type wgs \
    --outdir results_reannotated
```

| Param | Effect |
|---|---|
| `--skip_spliceai_squirls` | Skips both. Asserts `--vcf` already carries their INFO fields. Also drops the requirement for the four SQUIRLS/SpliceAI resource params, since nothing that runs reads them |
| `--precomputed_cadd` / `--precomputed_cadd_tbi` | A prescored CADD table (`#Chrom Pos Ref Alt RawScore PHRED`, bgzipped + `tabix -s 1 -b 2 -e 2`). Both or neither. On its own it *prescores*: variants found in it take its scores and only the remainder goes to `CADD.sh`, so `--cadd_data_dir` stays required. Together with `--skip_spliceai_squirls` (as above) it passes through untouched -- nothing is scored, and `--cadd_data_dir` is not required |

Prescoring accepts either contig naming -- the callset's (a previous run's `01_cadd/`
output, `chr1`/`chrM`) or CADD's own Ensembl naming (raw `CADD.sh` output, `1`/`MT`) --
and republishes the merged prescored + newly-scored table to `01_cadd/`, so each run's
table is a valid `--precomputed_cadd` for the next. The passthrough mode assumes the
table is complete **and** callset-named; neither is verified, both are warned about.

Three things worth knowing:

- **CADD is never skipped by a bare `--skip_cadd`.** ANNOVAR needs the table either way, so
  a flag with nothing behind it could not work: supply a table, and pair it with
  `--skip_spliceai_squirls` when nothing new should be scored.
- **`build_cadd_humandb` still runs.** Reformatting and indexing the table into ANNOVAR's
  `humandb/` layout is a different job from scoring a genome, and cheap next to it.
- **`--skip_spliceai_squirls` is not verified.** Nothing checks that the VCF really carries
  those INFO fields; a plain VCF produces empty splice columns rather than an error, because
  `add_splice_scores` reads absent INFO as absent. The run logs a warning naming the
  assumption. `03_annotated_vcf` is not republished on a skipped run -- you already have it.

`-resume` is not a substitute. It only works within one work directory with unchanged inputs,
and changing `--clinvar_date` for a new ClinVar release is exactly the case that invalidates
the ANNOVAR tasks while everything upstream should be reused.

Covered by [`../../tests/test_skip_params.py`](../../tests/test_skip_params.py), which runs
the pipeline against fake tools and asserts on which processes Nextflow actually submitted.

## Notes

- Every process emits a `versions.yml`; `mergeVersions` collects and deduplicates them, so a
  process that ran 24 per-interval tasks appears once.
- Memory is declared as `N.GB * task.attempt` throughout, so a too-low estimate self-corrects
  on retry instead of failing the run. The per-process numbers were split out of the original
  monolith's flat 50 GB and are estimates -- re-tune from the first real run's trace.
- `build_cadd_humandb` runs once per run and broadcasts, not once per interval.
- The annotation modules under `modules/annotations/` come from an upstream pipeline,
  with each file's header naming its local divergence. `annovar.nf` is a full rewrite of
  `wf`'s single monolithic per-interval process into six; `cadd.nf`, `spliceai.nf` and
  `squirls.nf` each carry a small local hunk.
- `bin/add_scores.py` is on PATH inside every task automatically -- Nextflow puts the entry
  script's `bin/` there. It is bundled, not configured, so there is no param for it. It
  deliberately does **not** exit on a multi-allelic record the way the upstream copy
  does: that pipeline decomposes in `merge_variants` upstream of the step, this one decomposes
  in `normalize_vcf` downstream of it, so multi-allelic records reach it on a normal run. See
  the comment in the record loop before syncing that hunk either way.
- `-stub-run` does not work here, and `-preview` hangs. Use `nextflow lint` and
  `nextflow inspect` -- but note both catch only parse and call-arity errors, not bad emit
  names or record shape mismatches.
  The tests fill that gap by running the real pipeline against fake tools on PATH:
  `./tests/run.sh` runs the whole suite, `./tests/run.sh skip_params` one test.

## Known gaps

The two that matter before trusting output:

1. No run has yet crossed the `squirls -> annovar` and `cadd -> build_cadd_humandb` record
   boundaries; declarations match by inspection only.
2. Single-interval output has not been byte-compared against the unmodified `wf`
   `run_annovar` it was rewritten from. That diff is the highest-value check available and
   needs compute2, containers and real data.
