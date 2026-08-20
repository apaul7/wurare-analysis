#!/usr/bin/env nextflow
nextflow.enable.types = true

/*
 * Merge, standardize, annotate and filter structural variants across an arbitrary mix of
 * SV VCFs -- any number of callers, any number of samples per VCF.
 *
 * Unlike annotate_snps, which is handed one already-merged VCF and justifies having no
 * samplesheet, this pipeline has to *build* the cohort VCF from heterogeneous inputs. That
 * needs a sheet, and it needs two: a multi-sample VCF cannot name one alignment, so
 * alignments live in their own sample-keyed sheet rather than overloading one file with
 * two record types.
 *
 * Phase 6: sheets parsed, preflight reported, inputs standardized, merged along both axes,
 * depth-annotated per sample, given gnomAD frequency and an AnnotSV report, then soft-tagged
 * in FILTER, with an optional Talos-schema tail. Nothing is ever removed.
 * The SVDB behaviours the design rests on are pinned by tests/test_svdb_assumptions.py.
 */

include { mergeVersions } from '../shared/modules/versions.nf'
include { to_sv_input; to_alignment; check_unique_rows } from './modules/other/samplesheet.nf'
include { prepare_svs } from './subworkflows/prepare_svs.nf'
include { merge_svs } from './subworkflows/merge_svs.nf'
include { depth_svs } from './subworkflows/depth_svs.nf'
include { qc_somalier } from './subworkflows/qc_somalier.nf'
include { annotate_cohort } from './subworkflows/annotate_cohort.nf'
include { tag_filters } from './modules/filter/tag_filters.nf'
include { talos_tail } from './subworkflows/talos_tail.nf'

workflow {
    main:
    if (!params.vcfs) error "params.vcfs (samplesheet of input SV VCFs) is required"

    // Required, not optional. It carries family structure, which
    // nothing else supplies, and its sex column is the cross-check that somalier's ploidy
    // call is measured against. Every cohort sample must have a row; preflight in prepare_svs
    // reports the ones that do not.
    if (!params.ped) {
        error "params.ped is required -- a 6-column PED (family, sample, father, mother, " +
              "sex, phenotype) supplying sample sex and family structure"
    }

    // Required only once alignments are in play, so that a run without them is not blocked
    // on a reference it will never read. CRAM cannot be decoded without the exact FASTA it
    // was written against, and the wrong one gives wrong depth rather than an error.
    if (params.alignments && !params.alignment_reference) {
        error "params.alignment_reference is required when params.alignments is set -- " +
              "CRAM decode needs the exact FASTA the alignments were written against"
    }
    if (params.alignments && !params.alignment_reference_index) {
        error "params.alignment_reference_index is required when params.alignments is set"
    }

    // Ancestry is all-or-none, like the annotation and Talos groups: half of it produces a
    // report missing the reference cohort it is predicted against, which is worse than none.
    ancestry_params = ['somalier_labels', 'somalier_1kg_dir']
    ancestry_configured = ancestry_params.findAll { name -> params[name] }
    if (ancestry_configured && ancestry_configured.size() != ancestry_params.size()) {
        error "somalier ancestry is partly configured: set both of " +
              "${ancestry_params.join(', ')} or neither. Missing: " +
              "${(ancestry_params - ancestry_configured).join(', ')}"
    }
    if (ancestry_configured && !params.somalier_sites) {
        error "params.somalier_sites is required for somalier ancestry -- ancestry is " +
              "predicted from the same extracted sites as relatedness"
    }

    // Both sheets are parsed eagerly, before any channel exists, for two reasons. A bad
    // row then fails at startup with the offending column named, the same way a missing
    // param does -- where parsing inside .map{} surfaces as an unhelpful
    // "Unexpected error [InvocationTargetException]", because error() raised in a closure
    // is wrapped before it reaches the console. And a sheet is small, so there is nothing
    // to stream.
    sv_inputs = channel.fromList(
        check_unique_rows(
            file(params.vcfs, checkIfExists: true)
                .splitCsv(header: true)
                .collect { row -> to_sv_input(row) }
        )
    )

    // Optional: a cohort with no alignments simply skips the depth stages, and preflight
    // lists the samples it will skip rather than dropping them quietly.
    alignments = params.alignments
        ? channel.fromList(
            file(params.alignments, checkIfExists: true)
                .splitCsv(header: true)
                .collect { row -> to_alignment(row) }
          )
        : channel.empty()

    // Shipped with the pipeline, not a param: it is source, not configuration.
    normalize_awk = file("${projectDir}/assets/normalize_records.awk", checkIfExists: true)
    stamp_awk = file("${projectDir}/assets/stamp_records.awk", checkIfExists: true)
    promote_awk = file("${projectDir}/assets/promote_caller_support.awk", checkIfExists: true)
    extract_depth_awk = file("${projectDir}/assets/extract_depth.awk", checkIfExists: true)
    merge_depth_awk = file("${projectDir}/assets/merge_depth.awk", checkIfExists: true)
    coverage_awk = file("${projectDir}/assets/check_annotsv_coverage.awk", checkIfExists: true)
    blob_awk = file("${projectDir}/assets/list_blob_info_keys.awk", checkIfExists: true)
    shard_awk = file("${projectDir}/assets/shard_vcf.awk", checkIfExists: true)
    concat_awk = file("${projectDir}/assets/concat_annotsv_tsv.awk", checkIfExists: true)
    filter_awk = file("${projectDir}/assets/tag_filters.awk", checkIfExists: true)
    check_somalier_awk = file("${projectDir}/../shared/assets/check_somalier.awk", checkIfExists: true)
    // Stands in for the ploidy table when somalier does not run, so tag_filters' Path input
    // is always a real file rather than a conditional channel shape. It names no sample, and
    // the awk then leaves chrX/chrY exempt -- the behaviour that predates ploidy awareness.
    ploidy_none = file("${projectDir}/assets/ploidy_none.tsv", checkIfExists: true)
    tsv_filter_awk = file("${projectDir}/assets/filter_annotsv_tsv.awk", checkIfExists: true)
    symbolic_alt_awk = file("${projectDir}/assets/symbolic_alt.awk", checkIfExists: true)
    talos_schema_awk = file("${projectDir}/assets/talos_schema.awk", checkIfExists: true)
    talos_keep_awk = file("${projectDir}/assets/talos_keep_info.awk", checkIfExists: true)
    talos_check_awk = file("${projectDir}/assets/check_talos_fields.awk", checkIfExists: true)

    ped = file(params.ped, checkIfExists: true)

    prepared = prepare_svs(sv_inputs, alignments, ped, normalize_awk, stamp_awk,
                           params.min_sv_size as Integer)

    // Axis A is cross-caller over the same samples and needs looser thresholds than axis
    // B's cross-sample assembly. Never one value for both just because the flag has one
    // name -- and when the cohort mixes joint and single-sample inputs, axis B takes the
    // looser cross-caller value too.
    merged = merge_svs(
        prepared.standardized, promote_awk,
        "${params.overlap_axis_a}", "${params.overlap_axis_b}", "${params.bnd_distance}",
        "${params.ins_distance}"
    )

    // Depth runs only where an alignment was supplied; with none, duphold never fires and
    // the cohort passes through untouched. Invoked once and its results bound to
    // locals -- a workflow, like a process, cannot be called twice in one context, and
    // reaching for .out after a call inside a ternary is the same mistake wearing a
    // different hat. The reference index is its own param rather than derived from the
    // FASTA: a derived index desyncs the moment someone overrides it.

    // Intersected against the cohort's own samples before it reaches duphold -- see
    // prepare_svs.nf. Preflight names the ignored rows.
    in_cohort = alignments.combine(prepared.cohort_samples)
        .filter { a, samples -> (samples as List).contains(a.sample) }
        .map { a, _samples -> a }

    if (params.alignments) {
        depth = depth_svs(
            merged.cohort, in_cohort,
            file(params.alignment_reference, checkIfExists: true),
            file(params.alignment_reference_index, checkIfExists: true),
            extract_depth_awk, merge_depth_awk, blob_awk
        )
        with_depth = depth.cohort
        depth_versions = depth.versions
    }
    else {
        with_depth = merged.cohort
        depth_versions = channel.empty()
    }

    // Phase 3b: sample QC from the same alignments -- ploidy for the depth filter, plus the
    // relatedness the Talos inheritance tail is entitled to have checked. Gated on its own
    // param rather than folded into the depth stage: --alignments without --somalier_sites
    // is a supported, degraded run, and the sex chromosomes then stay exempt exactly as they
    // were. Preflight already names what it skipped; this adds one line for the same reason.
    if (params.alignments && params.somalier_sites) {
        qc = qc_somalier(
            in_cohort,
            merged.cohort.first().map { c -> c.label },
            ped,
            file(params.somalier_sites, checkIfExists: true),
            file(params.alignment_reference, checkIfExists: true),
            file(params.alignment_reference_index, checkIfExists: true),
            check_somalier_awk,
            // Both are read only when the ancestry stage runs; ploidy_none stands in for the
            // unset case because a Path input cannot be null.
            ancestry_configured ? file(params.somalier_labels, checkIfExists: true)
                                : ploidy_none,
            ancestry_configured ? file(params.somalier_1kg_dir, checkIfExists: true)
                                : ploidy_none,
            ancestry_configured as Boolean
        )
        ploidy_tsv = qc.ploidy
        // The PED the rest of the pipeline uses: the operator's family structure, with sex
        // taken from the data. The original file is never modified.
        effective_ped = qc.inferred_ped
        qc_versions = qc.versions
        qc_samples = qc.samples
        qc_pairs = qc.pairs
        qc_html = qc.html
        qc_ancestry = qc.ancestry
        qc_ancestry_cohort = qc.ancestry_cohort
        qc_ancestry_html = qc.ancestry_html
        qc_inferred_ped = qc.inferred_ped
        // Published only when it was measured. The placeholder that stands in for it
        // otherwise is an asset, not a result, and publishing it would read as a ploidy
        // table whose every sample came back undetermined.
        qc_ploidy_out = qc.ploidy
    }
    else {
        ploidy_tsv = channel.value(ploidy_none)
        effective_ped = channel.value(ped)
        qc_versions = channel.empty()
        qc_samples = channel.empty()
        qc_pairs = channel.empty()
        qc_html = channel.empty()
        qc_ancestry = channel.empty()
        qc_ancestry_cohort = channel.empty()
        qc_ancestry_html = channel.empty()
        qc_inferred_ped = channel.empty()
        qc_ploidy_out = channel.empty()
    }

    // Annotation is conditional as a whole, not guarded with a hard error. This is a
    // deliberate departure from annotate_snps' "required param" pattern, and the reason is
    // that the two pipelines differ in what they are for: annotate_snps annotates the VCF it
    // is handed, so a missing resource means it cannot do its job. Here Phases 1-3 already
    // produce a merged, depth-annotated cohort VCF that stands on its own, and the AnnotSV
    // bundle is multiple gigabytes of site-local install. Failing the whole run at startup
    // because that bundle is absent would make the useful half unreachable.
    //
    // All four are required together once any is set -- a half-configured annotation stage
    // is worse than none, because it silently produces a report missing a whole database.
    annotation_params = ['gnomad_sv_vcf', 'gnomad_cnv_vcf',
                         'annotsv_annotations_dir', 'knotannotsv_config']
    configured = annotation_params.findAll { name -> params[name] }

    if (configured && configured.size() != annotation_params.size()) {
        error "annotation is partly configured: set all of " +
              "${annotation_params.join(', ')} or none. Missing: " +
              "${(annotation_params - configured).join(', ')}"
    }

    if (configured) {
        annotated = annotate_cohort(
            with_depth,
            file(params.gnomad_sv_vcf, checkIfExists: true),
            file(params.gnomad_cnv_vcf, checkIfExists: true),
            file(params.annotsv_annotations_dir, checkIfExists: true),
            file(params.knotannotsv_config, checkIfExists: true),
            coverage_awk,
            blob_awk,
            shard_awk,
            concat_awk,
            "${params.annotsv_shard_bytes}",
            "${params.genome_build}",
            "${params.annotsv_bundle_version}",
            "${params.query_overlap}", "${params.query_bnd_distance}",
            "${params.annotsv_drop_info}",
            "${params.gnomad_sv_in_occ}", "${params.gnomad_sv_in_frq}",
            "${params.gnomad_cnv_in_occ}", "${params.gnomad_cnv_in_frq}",
            tsv_filter_awk,
            "${params.tsv_filter_rare_af}",
            "${params.tsv_filter_tier1_acmg}", "${params.tsv_filter_tier2_acmg}",
            "${params.tsv_filter_pli_min}", "${params.tsv_filter_loeuf_max}",
            "${params.tsv_filter_rank_min}",
            filter_awk,
            ploidy_tsv,
            "${params.filter_pop_af}", "${params.filter_internal_af}",
            "${params.filter_del_dhffc}", "${params.filter_dup_dhbfc}",
            "${params.filter_min_callers}"
        )
        annotated_out = annotated.vcf
        tagged_out = annotated.tagged
        // tag_filters ran inside annotate_cohort, so its versions are already in that
        // subworkflow's mix. Nothing to add here.
        tag_versions = channel.empty()
        annotsv_tsv_out = annotated.tsv
        annotsv_tier1_out = annotated.tsv_tier1
        annotsv_tier2_out = annotated.tsv_tier2
        annotsv_coverage_out = annotated.coverage
        annotsv_report_out = annotated.report
        annotation_versions = annotated.versions
    }
    else {
        annotated_out = channel.empty()
        annotsv_tsv_out = channel.empty()
        annotsv_tier1_out = channel.empty()
        annotsv_tier2_out = channel.empty()
        annotsv_coverage_out = channel.empty()
        annotsv_report_out = channel.empty()
        annotation_versions = channel.empty()
    }

    // Tagging runs inside annotate_cohort when annotation ran, so `tagged_out` is already
    // bound. Without gnomAD there is no annotate_cohort and no COMMON_GNOMAD to write, so the
    // depth VCF gets tagged here instead. Same process, two call sites, which is legal: a
    // process cannot be invoked twice in ONE workflow context, and these are two.
    if (!configured) {
        final_vcf = with_depth.first()
        tagged = tag_filters(
            final_vcf.map { c -> c.label },
            final_vcf.map { c -> c.vcf },
            final_vcf.map { c -> c.tbi },
            final_vcf.map { c -> c.joint },
            filter_awk,
            ploidy_tsv,
            "${params.genome_build}",
            "${params.filter_pop_af}", "${params.filter_internal_af}",
            "${params.filter_del_dhffc}", "${params.filter_dup_dhbfc}",
            "${params.filter_min_callers}"
        )
        tagged_out = tagged.out
        tag_versions = tagged.versions
    }

    // Phase 6, optional and separate. Talos is one consumer, not the goal -- every
    // other consumer is already served by 05_filter, so this converts the cohort VCF into
    // the only schema Talos accepts and skipping it costs nothing else.
    //
    // Both inputs required together: PREDICTED_LOF needs the GTF and the BED, and a partial
    // configuration would emit a VCF Talos reads and silently under-reports from, which is
    // worse than not producing one. The PED used to be the third member of this group; it is
    // now a pipeline-level requirement, so the tail simply uses it.
    talos_params = ['protein_coding_gtf', 'noncoding_bed']
    talos_configured = talos_params.findAll { name -> params[name] }

    if (talos_configured && talos_configured.size() != talos_params.size()) {
        error "the Talos tail is partly configured: set all of ${talos_params.join(', ')} " +
              "or none. Missing: ${(talos_params - talos_configured).join(', ')}"
    }

    if (talos_configured) {
        talos = talos_tail(
            tagged_out,
            file(params.protein_coding_gtf, checkIfExists: true),
            file(params.noncoding_bed, checkIfExists: true),
            // The inferred PED, not params.ped: sex-stratified AF is exactly the consumer
            // that must not read a sex field this pipeline elsewhere declares untrustworthy.
            // Identical to params.ped when somalier did not run.
            effective_ped,
            symbolic_alt_awk, talos_schema_awk, talos_keep_awk, talos_check_awk,
            "${params.gnomad_pop}"
        )
        talos_vcf_out = talos.vcf
        talos_report_out = talos.report
        talos_versions = talos.versions
    }
    else {
        talos_vcf_out = channel.empty()
        talos_report_out = channel.empty()
        talos_versions = channel.empty()
    }

    software_versions = mergeVersions(
        prepared.versions.mix(merged.versions, depth_versions, qc_versions,
                              annotation_versions, tag_versions, talos_versions))

    publish:
    standardized_vcfs = prepared.standardized.map { g -> g.entry.vcf }
    standardized_tbis = prepared.standardized.map { g -> g.entry.tbi }
    filter_counts     = prepared.filter_counts
    preflight_txt     = prepared.preflight
    cohort_vcf        = merged.cohort.map { c -> c.vcf }
    cohort_tbi        = merged.cohort.map { c -> c.tbi }
    depth_vcf         = with_depth.map { c -> c.vcf }
    depth_tbi         = with_depth.map { c -> c.tbi }
    annotated_vcf     = annotated_out.map { c -> c.vcf }
    annotated_tbi     = annotated_out.map { c -> c.tbi }
    annotsv_tsv       = annotsv_tsv_out
    annotsv_tier1     = annotsv_tier1_out
    annotsv_tier2     = annotsv_tier2_out
    annotsv_coverage  = annotsv_coverage_out
    annotsv_report    = annotsv_report_out
    tagged_vcf        = tagged_out.map { c -> c.vcf }
    tagged_tbi        = tagged_out.map { c -> c.tbi }
    talos_vcf         = talos_vcf_out.map { c -> c.vcf }
    talos_tbi         = talos_vcf_out.map { c -> c.tbi }
    talos_fields      = talos_report_out
    qc_ploidy         = qc_ploidy_out
    qc_inferred_ped_f = qc_inferred_ped
    qc_samples_tsv    = qc_samples
    qc_pairs_tsv      = qc_pairs
    qc_relate_html    = qc_html
    // The full table (cohort + the 1kg reference set it was predicted against) and the same
    // table cut down to this run's samples.
    qc_ancestry_tsv   = qc_ancestry
    qc_ancestry_cohort_tsv = qc_ancestry_cohort
    qc_ancestry_page  = qc_ancestry_html
    pipeline_versions = software_versions
}
output {
    standardized_vcfs { path "01_prepare" }
    standardized_tbis { path "01_prepare" }
    filter_counts     { path "01_prepare" }
    // The routing decisions that used to go to stdout: which inputs grouped with which, which
    // samples have no alignment, which have no PED row. At cohort scale this is pages long.
    preflight_txt     { path "01_prepare" }
    cohort_vcf        { path "02_merge" }
    cohort_tbi        { path "02_merge" }
    // With no alignments this is the same file as 02_merge's: the depth stage ran and did
    // nothing, which is a real outcome and better published than silently absent.
    depth_vcf         { path "03_depth" }
    depth_tbi         { path "03_depth" }
    annotated_vcf     { path "04_annotate" }
    annotated_tbi     { path "04_annotate" }
    annotsv_tsv       { path "04_annotate" }
    // Subsets of the file above, published beside it and never in place of it. tier1 is
    // routinely empty -- see assets/filter_annotsv_tsv.awk for why that is the expected
    // result rather than a failure.
    annotsv_tier1     { path "04_annotate" }
    annotsv_tier2     { path "04_annotate" }
    annotsv_coverage  { path "04_annotate" }
    annotsv_report    { path "04_annotate" }
    tagged_vcf        { path "05_filter" }
    tagged_tbi        { path "05_filter" }
    talos_vcf         { path "06_talos" }
    talos_tbi         { path "06_talos" }
    talos_fields      { path "06_talos" }
    // Sample QC. Published even though only ploidy.tsv feeds another stage: the relatedness
    // tables and the HTML are the answer to "does this PED describe this cohort", which
    // nothing downstream can ask on its own, and the ploidy table is what a reviewer needs to
    // interpret a chrX call the depth filter did or did not tag.
    qc_ploidy         { path "07_qc" }
    qc_inferred_ped_f { path "07_qc" }
    qc_samples_tsv    { path "07_qc" }
    qc_pairs_tsv      { path "07_qc" }
    qc_relate_html    { path "07_qc" }
    qc_ancestry_tsv   { path "07_qc" }
    qc_ancestry_cohort_tsv { path "07_qc" }
    qc_ancestry_page  { path "07_qc" }
    pipeline_versions { path "pipeline_info" }
}
