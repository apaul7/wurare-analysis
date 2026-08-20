# annotate_svs

Takes an arbitrary mix of SV VCFs — any number of callers, any number of samples per VCF —
merges them into one cohort callset, adds depth and population evidence, annotates, and
soft-tags for filtering. Sibling to [`annotate_snps`](../annotate_snps/README.md), which
annotates one VCF it is handed; this one has to *build* the VCF first, and that is where the
hard decisions live.

**Needs Nextflow 26.04.0 or newer**, for the same reason as `annotate_snps`: every module sets
`nextflow.enable.types = true`, a preview flag absent from 25.10.x, and the manifest's
`nextflowVersion` cannot warn you because Nextflow compiles the modules before it checks that
field.

```bash
nextflow run nf/annotate_svs/main.nf \
    -profile compute2,wurare \
    --vcfs vcfs.csv \
    --ped cohort.ped \
    --alignments alignments.csv \
    --alignment_reference /path/to/reference.fa \
    --alignment_reference_index /path/to/reference.fa.fai \
    --outdir results
```

`--alignments` requires both `--alignment_reference` and `--alignment_reference_index`; the
run is rejected at startup otherwise. Drop all three to skip the depth stages entirely --
that is a supported run, not a degraded one, and preflight lists the samples it skipped.

## Input

**Two sheets, not one.** A multi-sample VCF cannot name a single alignment, and folding
alignment-only rows into the VCF sheet makes one file carry two record types with mutually
exclusive required columns.

`--vcfs` — required:

```csv
sample_set,caller,joint,vcf,tbi
FAM01,manta,false,/path/FAM01.manta.vcf.gz,/path/FAM01.manta.vcf.gz.tbi
FAM01,delly,false,/path/FAM01.delly.vcf.gz,/path/FAM01.delly.vcf.gz.tbi
COHORT_A,dragen,true,/path/cohortA.vcf.gz,/path/cohortA.vcf.gz.tbi
```

| Column | Notes |
|---|---|
| `sample_set` | Provenance label. `[A-Za-z0-9._-]` only — it reaches an output filename. **Not** the merge grouping key; that is derived from VCF headers |
| `caller` | `[A-Za-z0-9._-]` only — it becomes the SVDB `--vcf file:tag` tag and is interpolated unquoted into the merge command line |
| `joint` | `true`/`false`. Cannot be inferred from sample count: a multi-sample VCF made by concatenating single-sample calls is not joint-called, and giving it joint priority would overwrite the genotypes the joint-caller priority rule exists to protect |
| `vcf` / `tbi` | Sample IDs are read from the VCF header, never restated here — a second source of truth drifts, and drift means one sample silently becoming two columns |

A `(sample_set, caller)` pair must be unique. Two such rows stage under the same filename and
Nextflow keeps one, silently — so a duplicate is rejected at startup.

`--alignments` — optional, keyed on sample ID:

```csv
sample,alignment,alignment_index
FAM01-01,/path/FAM01-01.cram,/path/FAM01-01.cram.crai
```

A sample absent from this sheet skips the depth stages and is **listed at preflight**, not
dropped quietly. A cohort with no alignments at all is a supported configuration. The reverse
is too: a row naming a sample this cohort does not contain is ignored and named at preflight,
so one site-wide `alignments.csv` reused across cohorts works.

`--ped` — **required**, standard 6-column PED, whitespace-separated:

```
FAM01	FAM01-01	FAM01-02	FAM01-03	1	2
FAM01	FAM01-02	0	0	1	1
FAM01	FAM01-03	0	0	2	1
```

Columns are family, sample, father, mother, sex (`1` male, `2` female, `0` unknown),
phenotype. It carries the family structure nothing else supplies, and its sex column is the
**cross-check** against the sex somalier measures from the alignments — where the two
disagree, the data wins and the disagreement is reported (see `07_qc/` below). Every cohort
sample should have a row; **preflight names the ones that do not** — they get no family
structure, and no PED sex to check the measured one against.

A cohort with no recorded pedigree is not stuck: run `annotate_snps` first with
`--alignments` and `--somalier_sites`, and pass its `<outdir>/04_qc/cohort.inferred.ped`
here — this pipeline's somalier stage then cross-checks it like any other PED.

## Params

| Param | Default | Notes |
|---|---|---|
| `--vcfs` | — | required |
| `--ped` | — | **required.** Sample sex and family structure; see above |
| `--alignments` | `null` | optional; enables Phase 3 |
| `--alignment_reference` | `null` | **required with `--alignments`.** CRAM cannot be decoded without the exact FASTA it was written against, and the wrong one gives wrong depth rather than an error |
| `--alignment_reference_index` | `null` | required with `--alignments`. Its own param, not derived as `<reference>.fai` — a derived index desyncs the moment someone overrides the FASTA |
| `--somalier_sites` | `null` | optional, used with `--alignments`. somalier sites VCF for GRCh38, this pipeline's only build: [sites.hg38.vcf.gz](https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz), listed with the [releases](https://github.com/brentp/somalier/releases). **This is what makes the depth filter ploidy-aware on chrX/chrY**, and it also produces the relatedness QC. Without it those chromosomes stay exempt from the depth tag |
| `--somalier_labels` / `--somalier_1kg_dir` | `null` | optional ancestry QC, **both or neither**, and both need `--somalier_sites`. Labels: [ancestry-labels-1kg.tsv](https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv); reference cohort: [1kg.somalier.tar.gz](https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz), unpacked — the param names the **directory**, not the tarball. Nothing filters on ancestry yet |
| `--min_sv_size` | `50` | size floor at standardization; BND exempt |
| `--overlap_axis_a` | `0.6` | cross-caller merge, same samples |
| `--overlap_axis_b` | `0.8` | cross-sample cohort assembly |
| `--bnd_distance` | `2000` | both merge axes |
| `--query_overlap` | `0.6` | `svdb query` — its own knob, not the merge one |
| `--query_bnd_distance` | `10000` | as above |
| `--gnomad_sv_vcf` | `null` | genome SV callset |
| `--gnomad_cnv_vcf` | `null` | exome-derived CNV callset — **a different file**; one `gnomad_sv` param is how the two get conflated, and the result is wrong AFs reported confidently |
| `--annotsv_annotations_dir` | `null` | AnnotSV bundle |
| `--annotsv_bundle_version` | `3.5` | The bundle's release, recorded in `versions.yml`. **Update by hand when the bundle changes** — nothing validates it. It cannot be derived: the unpacked tree records no version, and AnnotSV reports its own version, not the bundle's |
| `--knotannotsv_config` | `null` | knotAnnotSV config YAML |
| `--genome_build` | `GRCh38` | passed to AnnotSV |
| `--filter_pop_af` | `0.01` | soft tag threshold |
| `--filter_internal_af` | `0.03` | soft tag threshold; see the small-N caveat below |
| `--filter_del_dhffc` | `0.7` | DEL depth support |
| `--filter_dup_dhbfc` | `1.3` | DUP depth support |
| `--filter_min_callers` | `2` | corroboration tag |
| `--protein_coding_gtf` | `null` | Talos tail; `gatk SVAnnotate` → `PREDICTED_LOF`, Talos's hard gate. Must be plain and hierarchically ordered — see below |
| `--noncoding_bed` | `null` | Talos tail; ditto |
| `--gnomad_pop` | `gnomad_v4.1` | Talos reads `{gnomad_pop}_sv_AF` and `{gnomad_pop}_sv_SVID` |
| `--ins_distance` | `25` | third SVDB merge threshold, both axes. 2.12.0's default; it was 50 in 2.8.1, so it is set explicitly rather than inherited |
| `--gnomad_sv_in_occ` / `--gnomad_sv_in_frq` | `AC` / `AF` | INFO keys the gnomAD **SV** file holds its count and frequency under |
| `--gnomad_cnv_in_occ` / `--gnomad_cnv_in_frq` | `SC` / `SF` | same for the gnomAD **CNV** file — different names, not a typo |
| `--annotsv_shard_bytes` | `1000000000` | max bytes per AnnotSV shard. Guards Tcl 8's hard 2 GiB single-value cap, not a memory limit |
| `--annotsv_drop_info` | `""` | INFO keys stripped before AnnotSV only; the emitted VCF keeps them. Excel caps a cell at 32767 chars |
| `--tsv_filter_rare_af` | `0.01` | tier-1 benign-SV AF ceiling |
| `--tsv_filter_tier1_acmg` / `--tsv_filter_tier2_acmg` | `4` / `3` | ACMG class floor per tier |
| `--tsv_filter_pli_min` / `--tsv_filter_loeuf_max` | `0.9` / `2` | tier-1 gene-constraint gate; either satisfies it |
| `--tsv_filter_rank_min` | `0.9` | tier-2 AnnotSV ranking-score floor |
| `--outdir` | `results` | |

The four annotation params are **all-or-none**. A half-configured annotation stage silently
produces a report missing a whole database, which is worse than not running it.

**Verify the four `gnomad_*_in_occ`/`_in_frq` keys against the files in hand before the first
real run.** `svdb` reads these by name; get them wrong and it skips every database variant,
warns per variant, and **exits 0** with an unannotated VCF — so `COMMON_GNOMAD` never fires
and the population filter silently does nothing.

```bash
bcftools view -h <gnomad_file.vcf.gz> | grep '^##INFO'
```

`--protein_coding_gtf` and `--noncoding_bed` are their own all-or-none pair, and setting
neither simply skips the Talos tail. `--ped` is not part of that group — it is required
pipeline-wide.

## Output

```
results/
  01_prepare/    per-input standardized VCFs + filter_counts.tsv
  02_merge/      cohort VCF: CALLER_SUPP, NCALLER, ALGORITHMS, AC/AN/AF
  03_depth/      + duphold DHFFC/DHBFC/DHFC/DHBZ per sample
  04_annotate/   + gnomAD frequencies; AnnotSV TSV, coverage report, HTML
  05_filter/     + soft FILTER tags and INTERNAL_AF (see the small-N note below)
  06_talos/      Talos-schema VCF + field report (optional tail)
  07_qc/         somalier: ploidy.tsv, cohort.inferred.ped, relatedness tables + HTML,
                 ancestry (optional). Needs --alignments and --somalier_sites
  pipeline_info/ software_versions.yml
```

**Filtering is soft.** `FILTER` carries `COMMON_GNOMAD`, `COMMON_INTERNAL`,
`DEPTH_UNSUPPORTED`, `NO_CALLER_SUPPORT`; no record is ever removed, and the process asserts
that. One run then serves Talos, manual review and a QC report without re-running the
expensive half.

## Preparing external inputs, and what is undocumented elsewhere

**somalier's sites file must match the alignment reference.** This pipeline is GRCh38 only, so
the file to download is
[sites.hg38.vcf.gz](https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz), listed
with the [releases](https://github.com/brentp/somalier/releases).
A sites file for a different build genotypes nothing; `check_somalier.awk` turns that into a hard
error naming the likely cause, rather than a cohort of confidently undetermined samples. No
indexing or preparation is needed — pass the `.vcf.gz` as it downloads.

Ancestry needs two more files, and only together:
[`ancestry-labels-1kg.tsv`](https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv)
and the unpacked
[`1kg.somalier.tar.gz`](https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz).
`--somalier_1kg_dir` names the **directory** the tarball unpacks to, not the tarball.

**`06_talos/` also loses the ALT sequence, for the same kind of reason.** `gatk SVAnnotate`
accepts only symbolic ALTs and breakends; a sequence-resolved one throws

```
java.lang.IllegalArgumentException: Unexpected ALT allele: TTTTTTCTTTCTTT...
Expected breakpoint or symbolic ALT allele representing a structural variant record.
```

and it throws rather than skipping, so a single such record ends the whole tail. Manta writes
these routinely — an insertion whose sequence it assembled, a deletion spelled out in full —
so it is ordinary caller output, and it cannot appear in any cohort built from symbolic-only
test data. `assets/symbolic_alt.awk` rewrites ALT to `<SVTYPE>` and cuts REF to its anchoring
base, on the SVAnnotate branch only. `04_annotate/` and `05_filter/` keep the sequence, which
is the more informative representation and which every other tool reads happily. A literal ALT
with no `SVTYPE` cannot be converted — nothing names the symbol — so it passes through and is
counted to stderr rather than guessed at from REF/ALT lengths.

**`06_talos/` is the one output where FILTER is not what the pipeline decided.** Talos runs

```python
mt.filter_rows(hl.is_missing(mt.filters) | (mt.filters.length() == 0))
```

*before* any category logic, so the pipeline's soft filter tags are not soft there — they delete the record.
Measured: a 2-record cohort tagged `COMMON_INTERNAL` loaded as **0 rows**, and
`NO_CALLER_SUPPORT` fires on any single-caller cohort, so a correct and fully annotated
callset could reach Talos and vanish with nothing reporting it.

So `assets/talos_schema.awk` moves the tags rather than obeying them: a record with a
non-empty FILTER gets `INFO/SOFT_FILTERS=<tags>` and `FILTER=PASS`, with semicolons turned
into commas because a semicolon separates INFO fields. A `.` FILTER is left alone — it means
no filter was applied, and hail already reads that as missing.

**Read `SOFT_FILTERS`, not FILTER, on this file only.** Every other published VCF — including
`05_filter/` — keeps its tags in FILTER where they belong. A `PASS` in `06_talos/` means "not
hidden from Talos", never "passed every threshold".

**Producing a Talos-ready VCF, verified end to end.** The pipeline has been run whole and its
`06_talos/cohort.talos.vcf.gz` loads in Talos 11.1.0. The run needs no cluster and no
annotation bundle — `--vcfs` and `--ped`, plus the tail's two resources:

```bash
nextflow run nf/annotate_svs/main.nf -c docker.config --vcfs vcfs.csv --ped cohort.ped --protein_coding_gtf resources/gatk_sv/gencode.v47.basic.protein_coding.canonical.gtf --noncoding_bed resources/gatk_sv/noncoding.sort.hg38.bed --outdir results
```

where `docker.config` is `docker { enabled = true }` — selecting no `-profile` gives a local
run, since everything cluster-specific lives inside `compute2` (see `shared/conf/compute2.config`).
On an arm64 Mac add `runOptions = "--platform linux/amd64"`; the biocontainers are amd64 only.

Two things about that output worth knowing before trusting it. **Skipping the annotation
params means no population frequency**: `gnomad_v4.1_sv_AF` is declared but valued on no
record, and Talos's AF filter is `hl.or_else(af, MISSING) < threshold`, so every variant
passes the population filter rather than being held to it. Set the four annotation params for
a real run. And **`06_talos/` publishes symlinks into `work/`** like every other stage, so
`docker run -v` on the published path fails to resolve them — dereference with `cp -L` first,
or mount the work directory too.

**Where the Talos tail's two resources come from.** `gatk SVAnnotate` produces `PREDICTED_LOF`,
which Talos hard gates on, and it needs a protein-coding GTF and a noncoding BED. Take both
from GATK-SV's own published hg38 resources rather than substituting a GENCODE release: the
GTF has to satisfy GATK's codec, not merely name the right genes (see below). These are the
paths `inputs/values/resources_hg38.json` in `broadinstitute/gatk-sv` names, and they live in
two *different* public buckets:

```bash
curl -sSfL -O https://storage.googleapis.com/gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/gencode.v47.basic.protein_coding.canonical.gtf
```
```bash
curl -sSfL -O https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed
```

| File | Size | sha256 (first 16) | Param |
|---|---|---|---|
| `gencode.v47.basic.protein_coding.canonical.gtf` | 289 MB | `f757171cabc01916` | `--protein_coding_gtf` |
| `noncoding.sort.hg38.bed` | 61 MB | `c05113c54fd22074` | `--noncoding_bed` |

Both are hg38/GRCh38, matching `params.genome_build`. They are site-local inputs like the
AnnotSV bundle — not committed, and `resources/` is gitignored.

**The GTF must be plain and hierarchically ordered, and two plausible files are not.** GATK's
GENCODE codec rejects a gzipped GTF outright (`no suitable codecs found` — decompress it), and
it requires each gene's records to run gene, then transcript, then exon. The neighbouring
`gencode.canonical_pc.gtf.gz` in the `gcp-public-data--broad-references` bucket is *coordinate*
sorted, so a UTR line precedes its own transcript and SVAnnotate dies partway through with a
bare `NullPointerException: ... because "transcript" is null` that names neither the file nor
the problem. Verified working: `PREDICTED_LOF=OR4F5` on a `chr1:65500-71000` deletion.

**`svtk` is not used.** It cannot run — its bundled `no_contigs_template.vcf` has a `#CHROM`
line with no sample column, which modern pysam rejects as invalid VCF, and it is loaded
unconditionally in every current biocontainer build. Normalization is
`assets/normalize_records.awk` instead, applied uniformly to every caller: `SVLEN` sign by
type, `END` resolved when absent or collapsed to `POS`, `SVTYPE` from a symbolic ALT.

**DUP-encoded-as-INS is not reconciled.** Collapsing the two is a claim about biology, not a
format repair — the caller may genuinely mean an insertion. Jasmine's `--dup_to_ins` exists
if it is ever wanted.

**Depth tags on chrX/chrY need `--somalier_sites`; without it those chromosomes stay exempt.**
The problem is the DUP path: `DHBFC` is normalized against GC-matched bins genome-wide, so a
hemizygous male chrX sits near 0.5 and a real 2-copy duplication there reads about 1.0 — under
the flat 1.3 threshold, so every true male chrX DUP would be tagged `DEPTH_UNSUPPORTED`. The
DEL path is unaffected: `DHFFC` is normalized against the variant's own flanks, which are
equally hemizygous, so it reads about 1.0 rather than 0.5.

With somalier, the DUP threshold is scaled by the copies the sample is expected to carry at
that locus (1.3 → 0.65 hemizygous) and the DEL threshold is left alone. PAR gets the diploid
thresholds for everyone; the GRCh38 coordinates ship in `assets/tag_filters.awk`, and any
other `--genome_build` falls back to exempting the sex chromosomes rather than guessing. A sample whose karyotype came back undetermined keeps the exemption too.

**Sex comes from the alignments, not from the PED.** PED sex cannot be assumed correct — a
mislabelled row and a swapped sample are indistinguishable from inside the pipeline — so
`07_qc/ploidy.tsv` records both, the data wins, and a disagreement is written into the VCF as
a `##SAMPLE_SEX=<...,agreement=DISAGREES>` header line as well as warned about on stderr.
`07_qc/cohort.inferred.ped` is the PED with column 5 replaced by the inferred sex; the Talos
tail's sex-stratified AF reads that rather than `--ped`, whose family structure is still the
operator's to state. `samtools idxstats` was evaluated for this and rejected.

**Internal AF has a hard floor at small N, and the tag is withheld there.** The smallest
non-zero frequency a cohort of N can express is `1/(2N)`. Below the size where that already
clears `--filter_internal_af`, `COMMON_INTERNAL` would mean "is a carrier" rather than "is
common", so it is not written at all and the run says so on stderr. `INTERNAL_AF` is still
written, so the number is available — only the judgement is withheld. Lean on the population
AF while the cohort is small.

`INTERNAL_AF` is **not** `INFO/AF`. `bcftools +fill-tags` computes `AF` over *called* alleles
only, and SVDB leaves a non-calling sample as `./.`, so a private variant scores `AF=0.5` at
any cohort size. `INTERNAL_AF` counts every sample, treating a no-call as reference.

**Genotypes are not joint.** duphold annotates depth; it does not regenotype. A `./.` in a
merged row means "not called by this input", and depth is what makes it interpretable — not
a genotype. Paragraph and GraphTyper2 would change that and are deliberately out of v1.

**Merge provenance is two different quantities.** `CALLER_SUPP`/`NCALLER` name the callers.
The cohort VCF's `SUPP_VEC`/`svdb_origin` name the *sample sets* that contributed, and its
bits are ordered alphabetically by tag rather than by input order — never read it
positionally.

**In-house frequency is not in v1.** `svdb build` needs a cohort to build from, and the
pipeline that produces that cohort is this one.

**BNDs are kept, not dropped — and their mate links are not maintained.** A breakend is
exempt from the size floor, from `END`/`SVLEN` normalization and from the depth tag, so
nothing about it is silently corrupted. But no merger re-pairs mates: after `svdb merge` a
breakend's `MATEID` may name a record that no longer exists under that ID. Translocations
and Manta's BND-encoded inversions therefore reach the output as individual breakends, and
consequence annotation on them is weak. Keeping them is the lesser cost — dropping them
loses every translocation — but do not read a BND row as a resolved event.

**`--priority` ordering.** Axis A takes `manta, smoove, cnvkit, cnvnator`, then any caller
not in that list alphabetically; axis B takes joint inputs first, ties broken by label. This
decides whose breakpoints and whose genotype survive a merged cluster, not merely the order
of arguments (verified by experiment during design). The axis-A order is site-local
knowledge, hardcoded in
`subworkflows/merge_svs.nf` — making it a param is an open item. Pinned by
`tests/test_priority_e2e.py`.

## Tests

```bash
conda env create -f ../../environment.yml
```

```bash
conda activate annotate-svs && ../../tests/run.sh
```

`run.sh` reports each test as ok/FAIL/skip and exits with the number of failures. A test that
cannot run — no Docker, no conda env — prints `SKIP:` and is counted separately, so an
all-skips run is not mistaken for a green one. Pass a name to run just one:
`../../tests/run.sh depth_merge`.

Most suites run the logic directly with no Nextflow and no container. Several do not:
`test_grouping_e2e.py`, `test_priority_e2e.py`, `test_depth_e2e.py` and `test_ploidy_e2e.py`
run the real pipeline under Docker, because what they check — the axis-A grouping key, the
`--priority` order at both axes, sample-to-alignment pairing, and the somalier ploidy table
reaching the filter — lives in Nextflow channel wiring and is not reproducible any other way.
`test_somalier_assumptions.py` needs Docker but not Nextflow: it runs somalier itself over a
toy genome to pin the columns the ploidy call reads. All of them skip cleanly when `nextflow`
or Docker is unavailable, so a green run without them is not a green run.

The rest exist because
every failure mode in this pipeline is silent: a wrong `SVLEN` sign, depth on the wrong
record, or a stale `ALGORITHMS` does not error, it just quietly reports something untrue.
`tests/test_svdb_assumptions.py` additionally pins SVDB's own behaviour — run it after any
SVDB upgrade, before the upgrade ships.

Two need a container rather than the conda env, and `conda activate annotate-svs` will not
get you either:

```bash
docker build --platform linux/amd64 -t wurare-svdb-tests tests/ && docker run --rm --platform linux/amd64 -v "$PWD:/w" -w /w wurare-svdb-tests python3 tests/test_svdb_assumptions.py
```

`test_svdb_assumptions.py` needs svdb *and* bcftools on one PATH, and bioconda ships svdb for
linux-64 only — so on an arm64 Mac `conda env create` cannot produce a working svdb, and the
suite silently skips itself on `shutil.which`. `tests/Dockerfile` supplies both.

```bash
TALOS_IMAGE=talos:11.1.0 python3 tests/test_talos_load.py
```

`test_talos_load.py` puts the schema question to Talos itself — it imports the real
`rearrange_annotations` and runs it over the tail's output, which is the only check here that
can catch Talos changing what it reads. Talos publish no image, so build theirs first
(`docker build -t talos:11.1.0 .` in their repo) and point `TALOS_IMAGE` at it. It skips, with
the reason, when the image is absent.

## Status

All six build phases are implemented.

Not yet verified against real data, all for want of site-local resources rather than design:
`svdb query` against a real gnomAD file, AnnotSV/knotAnnotSV (annotation bundle), and the
Talos tail (`gatk SVAnnotate` needs a protein-coding GTF and noncoding BED). The thresholds
— `--overlap_axis_a`/`_b` and the filter defaults — are starting points, not measurements;
tune them against a real cohort.
