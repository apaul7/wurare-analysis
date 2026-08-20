nextflow.enable.types = true
include { AlignmentEntry } from '../modules/other/samplesheet.nf'
include { somalier_extract; somalier_relate; somalier_ancestry;
          check_somalier } from '../../shared/modules/qc/somalier.nf'

// Phase 3b: sample QC from the alignments -- ploidy, relatedness, optionally ancestry.
//
// Deliberately NOT part of depth_svs, even though both consume the same alignments channel.
// The relatedness half is worth having on its own, the ploidy table is consumed by Phase 5
// rather than Phase 3, and duphold is the expensive stage this must not be entangled with.
workflow qc_somalier {
    take:
    alignments: Channel<AlignmentEntry>
    label: Channel<String>
    ped: Path
    sites: Path
    reference: Path
    reference_index: Path
    check_awk: Path
    // Ancestry is optional and all-or-none: the labelled reference cohort's directory and its
    // labels file. main.nf decides whether the stage runs at all.
    ancestry_labels: Path
    ancestry_dir: Path
    run_ancestry: Boolean

    main:
    // Records are destructured before a process sees them: one in an `input:` block has no
    // value-based hash, so its task hash changes every run and -resume never matches.
    parts = alignments.multiMap { a ->
        sample: a.sample
        alignment: a.alignment
        index: a.index
    }
    extracted = somalier_extract(parts.sample, parts.alignment, parts.index,
                                 sites, reference, reference_index)

    // Sorted before collecting, for the reason depth_svs and annotate_cohort also sort:
    // Nextflow emits collected items in COMPLETION order, so an unsorted collect changes the
    // task hash between runs and silently defeats -resume. The per-sample files are keyed on
    // sample name, so order cannot change the result -- only the hash.
    collected = extracted.out.toSortedList { a, b -> a.name <=> b.name }

    // [ped], not ped: somalier_relate takes the PED as an optional single-element list, because
    // annotate_snps can run without one. This pipeline always has one -- params.ped is required
    // in main.nf, and the inferred PED it produces is consumed downstream.
    // infer: false, always -- this pipeline requires an operator PED and only checks it, so
    // inference is never wanted here.
    related = somalier_relate(label, collected, [ped], false)
    checked = check_somalier(related.samples, related.pairs, ped, check_awk)

    if (run_ancestry) {
        ancestry = somalier_ancestry(label, collected, ancestry_labels, ancestry_dir)
        ancestry_out = ancestry.out
        ancestry_cohort_out = ancestry.cohort_out
        ancestry_page = ancestry.html
        ancestry_versions = ancestry.versions
    }
    else {
        ancestry_out = channel.empty()
        ancestry_cohort_out = channel.empty()
        ancestry_page = channel.empty()
        ancestry_versions = channel.empty()
    }

    // .first() on the per-sample versions: every extract task reports the same somalier
    // build, and N identical entries would just pad software_versions.yml.
    //
    // Collected in main: a chained .mix() inside an emit: expression does not resolve under
    // typed syntax ("No such variable").
    all_versions = extracted.versions.first()
        .mix(related.versions, checked.versions, ancestry_versions)

    emit:
    // Read by tag_filters at both of its call sites, so it has to be a value channel.
    ploidy: Channel<Path> = checked.ploidy.first()
    // The PED with column 5 replaced by the inferred sex. Family structure is untouched --
    // the operator states that, somalier only checks it.
    inferred_ped: Channel<Path> = checked.inferred_ped.first()
    samples: Channel<Path> = related.samples
    pairs: Channel<Path> = related.pairs
    html: Channel<Path> = related.html
    // The full table, this cohort's samples alongside the 1kg reference set they were
    // predicted against.
    ancestry: Channel<Path> = ancestry_out
    // The same table cut down to this run's samples -- what you actually read.
    ancestry_cohort: Channel<Path> = ancestry_cohort_out
    ancestry_html: Channel<Path> = ancestry_page
    versions: Channel<Path> = all_versions
}
