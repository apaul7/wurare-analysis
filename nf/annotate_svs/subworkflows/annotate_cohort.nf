nextflow.enable.types = true
include { MergedVcf } from '../modules/merge/svdb_merge.nf'
include { svdb_query as query_sv; svdb_query as query_cnv } from '../modules/annotations/annotate.nf'
include { compress_index as pack_sv; compress_index as pack_cnv } from '../modules/annotations/annotate.nf'
include { annotsv_shards; annotsv; concat_annotsv } from '../modules/annotations/annotate.nf'
include { check_coverage; knotannotsv } from '../modules/annotations/annotate.nf'
include { filter_tsv as filter_tier1; filter_tsv as filter_tier2 } from '../modules/annotations/annotate.nf'
include { tag_filters } from '../modules/filter/tag_filters.nf'

// Phase 4: population frequency onto the VCF, then AnnotSV as the terminal report.
//
// Named annotate_cohort rather than annotate_svs so it does not collide with the pipeline's
// own name.
//
// gnomAD SV (genome) and gnomAD CNV (exome-derived) are two separate files and two separate
// params. Conflating them produces wrong AFs, confidently reported -- which is why they get
// distinct output tags here rather than a shared "gnomad" one, and why there is no single
// `gnomad_sv` param anywhere in this pipeline.
workflow annotate_cohort {
    take:
    cohort: Channel<MergedVcf>
    gnomad_sv: Path
    gnomad_cnv: Path
    annotations_dir: Path
    knot_config: Path
    coverage_awk: Path
    blob_awk: Path
    shard_awk: Path
    concat_awk: Path
    // Max bytes per AnnotSV shard. See annotsv_shards: this exists because Tcl 8 caps a
    // single value at 2 GiB, not because of any memory limit.
    annotsv_shard_bytes: String
    genome_build: String
    // Recorded in versions.yml. Not derivable from the bundle -- see the annotsv process.
    annotsv_bundle_version: String
    overlap: String
    bnd_distance: String
    // The INFO keys each database holds its counts under. Two pairs and not one: the SV and
    // CNV releases do not use the same names, and getting them wrong is silent -- svdb skips
    // every database variant and still exits 0.
    // INFO keys dropped before AnnotSV only -- a reporting concern, see the process. SVDB's
    // per-input blobs are dropped there too, unconditionally and on that branch alone: the
    // VCF emitted below still carries them, so per-input provenance stays recoverable.
    annotsv_drop_info: String
    sv_in_occ: String
    sv_in_frq: String
    cnv_in_occ: String
    cnv_in_frq: String
    // Report-side filter thresholds. Every criterion is a knob so that "what did this file
    // actually filter on" is answerable from the file itself six months from now -- the awk
    // writes the active values into a header comment line on each filtered TSV.
    tsv_filter_awk: Path
    tsv_rare_af: String
    tsv_tier1_acmg: String
    tsv_tier2_acmg: String
    tsv_pli_min: String
    tsv_loeuf_max: String
    tsv_rank_min: String
    // Soft-filter tagging runs inside this subworkflow -- see the tag_filters call below.
    filter_awk: Path
    // Per-sample chrX/chrY copy numbers; assets/ploidy_none.tsv when somalier did not run.
    ploidy: Channel<Path>
    filter_pop_af: String
    filter_internal_af: String
    filter_del_dhffc: String
    filter_dup_dhbfc: String
    filter_min_callers: String

    main:
    // .first() makes this a value channel. `cohort` is read four times below, and a QUEUE
    // channel can only be consumed once -- the second reader gets nothing and the workflow
    // stalls or errors. This path had never run (it needs gnomAD and the AnnotSV bundle),
    // so nothing caught it.
    coh = cohort.first()
    label = coh.map { c -> c.label }
    joint = coh.map { c -> c.joint }
    // Records are destructured before they reach a process: one in an `input:` block is
    // hashed by object identity, so the task hash changes every run and -resume never
    // matches. They stay records in the plumbing here and in `emit:`.
    coh_vcf = coh.map { c -> c.vcf }
    coh_tbi = coh.map { c -> c.tbi }

    // Chained rather than parallel: each query adds its own INFO keys to the VCF the
    // previous one produced, so the carrier accumulates annotation instead of forking.
    // Aliased includes because a process cannot be invoked twice in one workflow.
    sv_queried = query_sv(label, coh_vcf, coh_tbi, gnomad_sv, sv_in_occ, sv_in_frq,
                          "gnomad_sv_OCC", "gnomad_sv_AF", overlap, bnd_distance)
    sv_vcf = pack_sv(label.map { l -> l + ".sv" }, sv_queried.out, joint)

    // .first() so the three reads below share one value channel. The label carried here is
    // "<label>.sv" from pack_sv above, and the CNV output filename is built from it.
    sv_packed = sv_vcf.out.first()
    cnv_queried = query_cnv(sv_packed.map { c -> c.label },
                            sv_packed.map { c -> c.vcf },
                            sv_packed.map { c -> c.tbi },
                            gnomad_cnv, cnv_in_occ, cnv_in_frq,
                            "gnomad_cnv_OCC", "gnomad_cnv_AF", overlap, bnd_distance)
    annotated = pack_cnv(label, cnv_queried.out, joint)
    // Read by annotsv AND check_coverage, so it has to be a value channel too.
    annotated_vcf = annotated.out.first()

    // Phase 5 tagging, here rather than after this subworkflow, because both ends of the
    // dependency chain live in here: COMMON_GNOMAD needs the AF the two queries above added,
    // and AnnotSV below must see the tags -- filter_annotsv_tsv.awk reads them out of the
    // TSV's FILTER column. Downstream, both tier filters tested tags that did not exist.
    //
    // Also normalizes FILTER: pass_size_filter admits -f PASS,. and the awk rewrites an
    // untagged "." to PASS, so tier 1 stops dropping those for a non-pathogenic reason.
    tagged = tag_filters(annotated_vcf.map { c -> c.label },
                         annotated_vcf.map { c -> c.vcf },
                         annotated_vcf.map { c -> c.tbi },
                         joint,
                         filter_awk, ploidy, genome_build,
                         filter_pop_af, filter_internal_af,
                         filter_del_dhffc, filter_dup_dhbfc, filter_min_callers)

    // Read by annotsv_shards AND check_coverage, so it has to be a value channel too.
    tagged_vcf = tagged.out.first()
    ann_label = tagged_vcf.map { c -> c.label }
    ann_vcf = tagged_vcf.map { c -> c.vcf }
    ann_tbi = tagged_vcf.map { c -> c.tbi }

    // AnnotSV last, because TSV is terminal: everything that wants a VCF is fed from
    // `annotated_vcf` above, not from here.
    //
    // Sharded, because AnnotSV is Tcl and Tcl 8 cannot hold a value over 2 GiB -- a
    // 300-sample cohort VCF is ~9.2 GB even after stripping. One task per shard, so the
    // scheduler decides the concurrency; no maxForks here.
    shards = annotsv_shards(ann_label, ann_vcf, ann_tbi, annotsv_drop_info, blob_awk,
                            shard_awk, annotsv_shard_bytes)
    per_shard = annotsv(shards.shards.flatten(), annotations_dir, genome_build,
                        annotsv_bundle_version)

    // Sorted before collecting, for the same reason depth_svs sorts its per-sample outputs:
    // Nextflow emits collected items in completion order, so an unsorted collect would put
    // the shards' rows in a different order on every run. The shard names are zero-padded and
    // ascending, so sorting on the name restores coordinate order.
    tsv = concat_annotsv(ann_label,
                         per_shard.out.toSortedList { a, b -> a.name <=> b.name },
                         concat_awk)
    // Read by check_coverage, knotannotsv and both filter tiers.
    tsv_value = tsv.out.first()

    coverage = check_coverage(ann_label, ann_vcf, ann_tbi, tsv_value, coverage_awk)
    report = knotannotsv(label, tsv_value, knot_config, genome_build)

    // Subsets of the report, not replacements for it: `tsv` above is emitted and published
    // unchanged. Aliased includes because a process cannot be invoked twice in one workflow.
    tier1 = filter_tier1(label, tsv_value, tsv_filter_awk, "1", tsv_rare_af,
                         tsv_tier1_acmg, tsv_pli_min, tsv_loeuf_max, tsv_rank_min)
    tier2 = filter_tier2(label, tsv_value, tsv_filter_awk, "2", tsv_rare_af,
                         tsv_tier2_acmg, tsv_pli_min, tsv_loeuf_max, tsv_rank_min)

    // .first() on the shard versions: every shard reports the same AnnotSV build, and N
    // identical entries would just pad software_versions.yml.
    all_versions = sv_queried.versions
        .mix(sv_vcf.versions, cnv_queried.versions, annotated.versions,
             tagged.versions,
             shards.versions, per_shard.versions.first(), tsv.versions,
             coverage.versions, report.versions,
             tier1.versions, tier2.versions)

    emit:
    // Two VCFs, deliberately distinct files. `vcf` is the gnomAD-annotated cohort before
    // tagging and is what 04_annotate publishes; `tagged` carries the FILTER tags and is what
    // 05_filter publishes and what AnnotSV actually read.
    vcf: Channel<MergedVcf> = annotated_vcf
    tagged: Channel<MergedVcf> = tagged_vcf
    tsv: Channel<Path> = tsv_value
    // Tier 1 is routinely EMPTY -- zero to a handful of rows per sample is the expected
    // outcome, not a failed run. Tier 2 is the working list, hundreds to thousands.
    tsv_tier1: Channel<Path> = tier1.out
    tsv_tier2: Channel<Path> = tier2.out
    coverage: Channel<Path> = coverage.out
    report: Channel<Path> = report.out
    versions: Channel<Path> = all_versions
}
