nextflow.enable.types = true
include { AlignmentEntry } from '../modules/other/samplesheet.nf'
include { MergedVcf } from '../modules/merge/svdb_merge.nf'
include { split_sample; duphold; recombine_depth } from '../modules/depth/duphold.nf'

// Phase 3: depth evidence per sample, recombined onto the cohort VCF.
//
// Runs only for samples that appear in the alignments sheet. A sample without one is not an
// error and is not dropped -- it comes out with "." for every depth field, and preflight
// has already listed it. A cohort with no alignments at all skips duphold entirely and
// passes the cohort VCF through unchanged, which is why recombine_depth tolerates an empty
// depth table rather than treating it as a failure.
workflow depth_svs {
    take:
    cohort: Channel<MergedVcf>
    alignments: Channel<AlignmentEntry>
    reference: Path
    reference_index: Path
    extract_awk: Path
    merge_awk: Path
    blob_awk: Path

    main:
    // The cohort VCF is one item and every duphold task needs it, so it is broadcast with
    // .first() -- a queue channel would be consumed by the first sample and the rest would
    // stall waiting for an item that never arrives.
    cohort_value = cohort.first()
    coh_label = cohort_value.map { c -> c.label }
    coh_vcf = cohort_value.map { c -> c.vcf }
    coh_tbi = cohort_value.map { c -> c.tbi }

    // Records are destructured before a process sees them: one in an `input:` block is
    // hashed by object identity, so its task hash changes every run and -resume never
    // matches it.
    //
    // Split, annotate, recombine. Three steps because the duphold image carries no
    // bcftools and duphold wants a single-sample VCF anyway.
    //
    // split_sample gets the sample name and nothing else about the alignment. Routing the
    // CRAM through it -- which is what carrying an AlignmentEntry on the sites record did --
    // cannot work: a Path rebuilt into an output record is a TaskPath pointing into that
    // task's work dir, and feeding one back into a process throws
    // UnsupportedOperationException out of the task hasher.
    sites = split_sample(coh_vcf, coh_tbi, alignments.map { a -> a.sample }, blob_awk)

    // So the alignment comes straight from the samplesheet channel and is JOINED to the
    // sites on the sample name. A keyed join, not a positional zip: duphold tasks finish out
    // of order, and pairing a sample with another sample's depth is silent and severe.
    // duphold re-checks the sites filename against the sample name on top of this.
    paired = sites.out
        .map { s -> tuple(s.sample, s.vcf, s.tbi) }
        .join(alignments.map { a -> tuple(a.sample, a.alignment, a.index) }, by: 0)

    parts = paired.multiMap { sample, vcf, tbi, alignment, index ->
        sample: sample
        vcf: vcf
        tbi: tbi
        alignment: alignment
        index: index
    }
    annotated_vcfs = duphold(parts.sample, parts.vcf, parts.tbi,
                             parts.alignment, parts.index,
                             reference, reference_index)

    // Sorted before collecting: the per-sample files are keyed on ID and sample name so
    // order cannot change the result, but an unsorted collect makes the task hash vary
    // between runs and silently defeats -resume. Same reasoning as the merge inputs in
    // merge_svs.nf.
    collected = annotated_vcfs.out.toSortedList { a, b -> a.name <=> b.name }

    annotated = recombine_depth(coh_label, coh_vcf, coh_tbi,
                                cohort_value.map { c -> c.joint },
                                collected, extract_awk, merge_awk)

    all_versions = sites.versions.mix(annotated_vcfs.versions, annotated.versions)

    emit:
    cohort: Channel<MergedVcf> = annotated.out
    versions: Channel<Path> = all_versions
}
