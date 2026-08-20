nextflow.enable.types = true
include { MergedVcf } from '../merge/svdb_merge.nf'

// Phase 4: population frequency, then consequence and ACMG.
//
// The shape of this stage is set by one fact: AnnotSV emits TSV, not VCF. That is a terminal
// format -- anything downstream that wants a VCF (Talos, IGV, a browser) has to be fed from
// a different branch. So the VCF stays the carrier through every VCF->VCF annotator and
// AnnotSV runs last. This is exactly what annotate_snps does with ANNOVAR; AnnotSV is this
// pipeline's ANNOVAR.

// gnomAD-SV frequency by reciprocal overlap. `svdb query` also serves the in-house frequency
// database later -- one tool, two databases -- but that DB needs a cohort to build from
// and the pipeline that produces it is this one, so it lands in v2.
//
// The thresholds here are NOT the merge thresholds. `svdb query` has its own --overlap and
// --bnd_distance with their own defaults, which is the third set in this pipeline after
// axis A and axis B. Reusing one value across all three because the flag has one name is the
// same mistake in a third place.
process svdb_query {
    container 'quay.io/biocontainers/svdb:2.12.0--py312hfcd9dac_0'
    cpus 1
    // The query database is loaded into memory and the gnomAD-SV genome VCF is large. This
    // is a sizing decision, not a default -- expect to tune it against the real file.
    memory { 16.GB * task.attempt }

    input:
    label: String
    vcf: Path
    tbi: Path
    db: Path
    // Which INFO keys to READ from the database, as against occ_tag/frq_tag, which name the
    // ones written onto the cohort VCF. svdb has no usable default: without these it skips
    // every database variant with a warning, exits 0 on an unannotated VCF, and then
    // COMMON_GNOMAD never fires and the population filter silently does nothing. On the
    // first real run it did not even get that far -- with nothing loaded, svdb crashed
    // indexing an empty array (IndexError, query_module.py:283). Per-file, because the SV
    // and CNV releases do not agree on the key names.
    in_occ: String
    in_frq: String
    occ_tag: String
    frq_tag: String
    overlap: String
    bnd_distance: String

    output:
    out: Path = file("${label}.${occ_tag}.vcf")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    svdb --query \\
        --query_vcf "${vcf.name}" \\
        --db "${db.name}" \\
        --in_occ "${in_occ}" \\
        --in_frq "${in_frq}" \\
        --out_occ "${occ_tag}" \\
        --out_frq "${frq_tag}" \\
        --overlap ${overlap} \\
        --bnd_distance ${bnd_distance} \\
        > "${label}.${occ_tag}.vcf"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    svdb: \$(svdb 2>&1 | head -n1 | sed 's/^usage: SVDB-//; s/:.*//')
    gnomad_db: ${db.name}
END_VERSIONS
    """
}

process compress_index {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    label: String
    vcf: Path
    joint: Boolean

    output:
    out: MergedVcf = new MergedVcf(
        label: label,
        joint: joint,
        vcf: file("${label}.annotated.vcf.gz"),
        tbi: file("${label}.annotated.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    bcftools view -Oz -o "${label}.annotated.vcf.gz" "${vcf.name}"
    bcftools index -t "${label}.annotated.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// The terminal report. Three things are deliberate here.
//
// The VCF this reads is a stripped copy, not the published one. AnnotSV puts INFO into the
// TSV verbatim, and a merged cohort's INFO is mostly SVDB's per-input blobs -- each one a
// whole source record as text, one per input, on every row. That is what makes the report
// unopenable in Excel rather than merely wide. They are dropped here and nowhere else: the
// cohort VCF published from 04_annotate still carries every blob, so the per-input
// provenance -- which input supplied which genotype -- is answerable from the VCF, which is
// where a provenance question belongs. Caller support survives the strip in its own right,
// because promote_caller_support.awk harvested CALLER_SUPP/NCALLER out of the blobs into
// typed INFO keys back at the merge.
//
// `-annotationMode both` is passed explicitly. It is also AnnotSV's default, but naming it
// here is the point: the mode decides what a row MEANS, and that has to be a stated choice
// rather than an inherited one. `both` emits full rows (one per SV) AND split rows (one per
// SV x gene) into one TSV, tagged in the Annotation_mode column.
//
// An earlier version passed `full`, on the reasoning that a reader expecting a per-SV table
// should not silently get a per-gene one. That was the wrong trade once the report gained
// filter tiers: the gene-level columns a rare-variant filter needs -- Location,
// Overlapped_CDS_percent, GnomAD_pLI, OMIM_morbid -- are per-transcript, and AnnotSV leaves
// them BLANK on full rows. Filtering on coding-ness against a full-only TSV is not strict,
// it is inert. So the published TSV carries both row types and every consumer states which
// it wants: the filter tiers keep split rows only, and check_coverage counts distinct IDs
// rather than rows precisely so that duplicate IDs across split rows cannot inflate it.
//
// The annotation bundle is a site-local param with no default. Its version is NOT captured:
// versions.yml records the AnnotSV binary only, so the bundle -- and the gnomAD-SV release it
// carries, which must not silently differ from the one svdb_query is pointed at -- goes
// unrecorded. Capturing both remains an open item.
// Strip once, then split, so AnnotSV never sees a value Tcl cannot hold.
//
// AnnotSV is Tcl, and Tcl 8 caps a SINGLE VALUE at 2^31-1 bytes. Even after the strip below,
// a 300-sample cohort VCF measured 9.24 GB -- ~300 sample columns across ~3M records -- and
// AnnotSV died in its first step ("VCF to BED") with "max size for a Tcl value (2147483647
// bytes) exceeded". That is a ceiling and not a sizing problem: no amount of memory moves it.
//
// `bcftools view -G` would have fixed it in two tokens by dropping the genotype columns, at
// the cost of "which sample carries this" disappearing from the report. Sharding keeps every
// column instead. Records go into fixed-size shards rather than per-contig ones because
// contigs are uneven -- chr1 alone is ~8% of the genome -- and an uneven shard has the same
// cliff behind it. The input is coordinate-sorted, so consecutive N-record blocks are
// contiguous coordinate ranges: contigs stay grouped and the assembled report stays in
// coordinate order anyway.
process annotsv_shards {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 8.GB * task.attempt }

    input:
    label: String
    vcf: Path
    tbi: Path
    // INFO keys to drop before AnnotSV, comma-separated, empty to keep everything. AnnotSV
    // copies INFO into the TSV verbatim, and smoove's PRPOS/PREND are long enough to push a
    // cell past Excel's 32767-character limit -- which makes the report unopenable for the
    // people it exists for. Dropped on this branch ONLY: the published cohort VCF keeps
    // every key, because the fix is a reporting concern and must not lose data.
    drop_info: String
    // Lists SVDB's per-input blobs, which are stripped here unconditionally and are the
    // larger half of the problem -- named keys alone did not make the report readable. Their
    // key names carry each input's tag, so they cannot be named in a param.
    blob_awk: Path
    shard_awk: Path
    // Max BYTES per shard, not records. Bytes-per-record scales with sample count -- ~12.6 KB
    // at 300 samples, nearly all genotype columns -- so a record count tuned for one cohort
    // is silently wrong for the next. A 200000-record default set from a guessed record count
    // produced a 2.52 GB shard on the first real cohort, 17% over the ceiling.
    shard_bytes: String

    output:
    shards: List<Path> = files('shard_*.vcf')
    versions: Path = file("versions.yml")

    script:
    def named = drop_info
        ? drop_info.split(',').collect { t -> "INFO/" + (t as String).trim() }.join(',')
        : ''
    """
    set -euo pipefail
    # The blobs are read off the header rather than listed, because their key names are the
    # input tags of this particular cohort. Empty on a cohort that never went through a merge.
    # Done ONCE here rather than per shard -- it is a full pass over the cohort VCF.
    blobs=\$(bcftools view -h "${vcf.name}" | awk -f "${blob_awk.name}")
    strip="${named}"
    if [ -n "\$blobs" ]; then strip="\${strip:+\$strip,}\$blobs"; fi

    if [ -n "\$strip" ]; then
        bcftools annotate -x "\$strip" -Oz -o "stripped.vcf.gz" "${vcf.name}"
    else
        cp "${vcf.name}" "stripped.vcf.gz"
    fi

    bcftools view -h "stripped.vcf.gz" > "header.vcf"
    n_shards=\$(bcftools view -H "stripped.vcf.gz" \\
        | awk -v MAXBYTES="${shard_bytes}" -v HDR="header.vcf" -f "${shard_awk.name}")

    # Every record must land in exactly one shard. A split that dropped records would produce
    # a report quietly missing variants, which check_annotsv_coverage.awk would only surface
    # much later and only if someone read it.
    before=\$(bcftools view -H "stripped.vcf.gz" | wc -l | tr -d ' ')
    after=\$(cat shard_*.vcf | grep -vc '^#' || true)
    if [ "\$before" != "\$after" ]; then
        echo "ERROR: sharding changed the record count: \$before -> \$after" >&2
        exit 1
    fi
    echo "sharded \$before records into \$n_shards shard(s)" >&2

    # The whole point of this process, so it is asserted rather than hoped for. A shard at or
    # over the Tcl ceiling would fail inside AnnotSV with a message that names neither the
    # shard nor the cause.
    for s in shard_*.vcf; do
        sz=\$(wc -c < "\$s" | tr -d ' ')
        if [ "\$sz" -ge 2147483647 ]; then
            echo "ERROR: \$s is \$sz bytes, at or over Tcl's 2147483647-byte value limit;" >&2
            echo "       AnnotSV would fail on it. Lower --annotsv_shard_bytes." >&2
            exit 1
        fi
    done

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// One task per shard. The output is named from the SHARD, not from a label carried alongside
// it, so there is no pairing to get wrong -- cf. the arrival-order defect that made
// promote_axis_a uncacheable and could mislabel a whole group.
process annotsv {
    container 'quay.io/biocontainers/annotsv:3.5.9--hdfd78af_0'
    cpus 2
    memory { 16.GB * task.attempt }

    input:
    shard: Path
    annotations_dir: Path
    genome_build: String
    // Operator-supplied, and NOT verified against the bundle -- nothing in the unpacked tree
    // records its release. AnnotSV reports its own version, never the bundle's, and the release
    // number exists only in the download tarball's name, which untar destroys. So this is a
    // label someone is responsible for keeping true, not a measurement.
    bundle_version: String

    output:
    out: Path = file("${shard.simpleName}.annotsv.tsv")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    AnnotSV \\
        -SVinputFile "${shard.name}" \\
        -annotationsDir "${annotations_dir.name}" \\
        -genomeBuild "${genome_build}" \\
        -annotationMode both \\
        -outputDir . \\
        -outputFile "${shard.simpleName}.annotsv.tsv"

    # AnnotSV writes nothing when a shard yields no annotatable records, which is legitimate.
    # concat_annotsv_tsv.awk tolerates an empty shard; a MISSING file would instead fail the
    # task on a missing output, so it is created rather than left absent.
    [ -f "${shard.simpleName}.annotsv.tsv" ] || : > "${shard.simpleName}.annotsv.tsv"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    AnnotSV: \$(AnnotSV -version 2>&1 | head -n1 | sed 's/^AnnotSV //')
    annotsv_bundle: ${bundle_version}
END_VERSIONS
    """
}

// Reassemble the shards into the single TSV the rest of the pipeline already consumes.
// The header handling is the whole job -- see assets/concat_annotsv_tsv.awk.
process concat_annotsv {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    label: String
    tsvs: List<Path>
    concat_awk: Path

    output:
    out: Path = file("${label}.annotsv.tsv")
    versions: Path = file("versions.yml")

    script:
    def shard_files = tsvs.collect { t -> t.name }.join(' ')
    """
    set -euo pipefail
    awk -f "${concat_awk.name}" ${shard_files} > "${label}.annotsv.tsv"

    # Rows in must equal rows out, headers aside. The awk drops every shard's header and
    # re-emits one; a bug there would silently lose a shard rather than error.
    in_rows=\$(cat ${shard_files} | grep -c . || true)
    n_hdrs=\$(for f in ${shard_files}; do [ -s "\$f" ] && echo x; done | wc -l | tr -d ' ')
    out_rows=\$(grep -c . "${label}.annotsv.tsv" || true)
    expected=\$(( in_rows - n_hdrs + 1 ))
    if [ "\$out_rows" != "\$expected" ]; then
        echo "ERROR: concatenation changed the row count: expected \$expected, got \$out_rows" >&2
        echo "       (\$in_rows non-empty lines in, \$n_hdrs shard(s) with a header)" >&2
        exit 1
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// Rare / potentially-damaging subsets of the report, tier 1 (strict) and tier 2 (loose).
// The unfiltered TSV is published beside these and is never replaced -- a view is allowed to
// drop rows only because the thing it is a view of is still there.
//
// bcftools' container purely for its awk, exactly as check_coverage does. Adding an image
// whose only job is to supply /usr/bin/awk buys nothing.
//
// The criteria live in assets/filter_annotsv_tsv.awk and every threshold is a param, for the
// same reason as everywhere else here: a stray dollar or quote in a triple-quoted Nextflow
// script shifts every later interpolation by one, silently.
process filter_tsv {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    label: String
    tsv: Path
    filter_awk: Path
    tier: String
    rare_af: String
    acmg_min: String
    pli_min: String
    loeuf_max: String
    rank_min: String

    output:
    out: Path = file("${label}.annotsv.tier${tier}.tsv")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    awk -v TIER="${tier}" -v RARE_AF="${rare_af}" -v ACMG_MIN="${acmg_min}" \\
        -v PLI_MIN="${pli_min}" -v LOEUF_MAX="${loeuf_max}" -v RANK_MIN="${rank_min}" \\
        -f "${filter_awk.name}" "${tsv.name}" \\
        > "${label}.annotsv.tier${tier}.tsv"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// Which records reached the TSV, and which did not -- by ID, never by row count.
process check_coverage {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 2.GB * task.attempt }

    input:
    label: String
    vcf: Path
    tbi: Path
    tsv: Path
    coverage_awk: Path

    output:
    out: Path = file("${label}.annotsv_coverage.txt")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    bcftools view "${vcf.name}" > "cohort.vcf"
    awk -f "${coverage_awk.name}" "cohort.vcf" "${tsv.name}" \\
        > "${label}.annotsv_coverage.txt"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// A cheap, genuinely useful last step for anyone reading results by hand: the AnnotSV TSV
// as an HTML report.
process knotannotsv {
    container 'quay.io/biocontainers/knotannotsv:1.1.5--hdfd78af_0'
    cpus 1
    memory { 8.GB * task.attempt }

    input:
    label: String
    tsv: Path
    config: Path
    genome_build: String

    output:
    out: Path = file("${label}.annotsv.html")
    versions: Path = file("versions.yml")

    script:
    // knotAnnotSV defaults to hg19 and takes UCSC names, while AnnotSV takes assembly names
    // -- the two have to be told the same thing in two vocabularies. Unpassed, the report
    // rendered GRCh38 coordinates with hg19 hyperlinks: a reviewer clicking a variant landed
    // on the same numbers in the wrong assembly, silently offset. The config YAML carries no
    // build key, so this flag is the only place it can be set.
    def ucsc_build = genome_build == "GRCh37" ? "hg19" : "hg38"
    """
    set -euo pipefail
    knotAnnotSV.pl \\
        --configFile "${config.name}" \\
        --genomeBuild "${ucsc_build}" \\
        --annotSVfile "${tsv.name}" \\
        --outDir . \\
        --outPrefix "${label}"

    if [ ! -f "${label}.annotsv.html" ]; then
        found=\$(find . -name '*.html' | head -n1)
        [ -n "\$found" ] && mv "\$found" "${label}.annotsv.html"
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    knotAnnotSV: \$(ls /usr/local/conda-meta | grep '^knotannotsv-' | head -n1 | cut -d- -f2)
END_VERSIONS
    """
}
