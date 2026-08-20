nextflow.enable.types = true
include { Reference } from '../../shared/modules/types.nf'
include { interval_tag; list_chroms; split_vcf } from '../../shared/modules/other/tools.nf'
include { cadd } from '../modules/annotations/cadd.nf'
include { spliceai } from '../modules/annotations/spliceai.nf'
include { squirls } from '../modules/annotations/squirls.nf'
include { build_cadd_humandb } from '../modules/annotations/build_cadd_humandb.nf'
include { annovar } from '../modules/annotations/annovar.nf'
include { merge_vcfs } from '../modules/other/merge_vcfs.nf'

// Split the input cohort VCF by interval, fan out to cadd/spliceai/squirls, then run the
// ANNOVAR chain per interval. Per-interval multianno.txt fragments are merged at the text
// level (header once + tail-concat) inside the annovar subworkflow -- the per-interval
// VCFs are never recombined.
//
// The split/fan-out half is ported from the upstream pipeline's annotate_variants.nf,
// which diverges after squirls: it goes on to vcfanno/VEP/variants_to_table, this goes to
// ANNOVAR. Comments on list_chroms, interval_tag and the .first() broadcast are carried
// over because they document real typed-preview constraints, not this pipeline's choices.
workflow annotate_snps {
    take:
    vcf_tbi: Channel<Path, Path>
    reference: Reference
    cadd_data_dir: Path
    squirls_config: Path
    squirls_jannovar: Path
    spliceai_precomputed_scores: Path
    spliceai_precomputed_tbi: Path
    cohort: String
    data_type: String
    clinvar_date: String
    // Already CLINVAR-substituted and preflighted against humandb/ in main.nf.
    annovar_protocols: String
    annovar_operations: String
    omim_xref: Path
    annovar_dir: Path
    splice_scores_script: Path
    // Empty list, or [cadd_tsv_gz, cadd_tbi]. A prescored CADD table -- e.g. a previous
    // run's 01_cadd output -- whose scores are copied so only uncovered variants are
    // re-scored; with skip_spliceai_squirls it is instead used as-is. Not a Tuple: typed
    // syntax has no null, and an empty list is how "absent" is expressed for an optional
    // file input.
    precomputed_cadd: List<Path>
    skip_spliceai_squirls: Boolean

    main:
    // Intervals are derived from what's actually in vcf_tbi (list_chroms), not a static
    // chr1..chrY list -- on targeted-panel data most chromosomes have zero variants, and
    // CADD/SpliceAI/SQUIRLS all error out on an empty VCF.
    chroms = list_chroms(vcf_tbi)
    chroms_out = chroms.out

    // Each interval is carried as (raw string, filename-safe tag) from here on. The tag is
    // computed once, at the point the interval enters the pipeline, so that every
    // downstream process names its outputs from `tag` while bcftools -r and SpliceAI's
    // region lookup keep using the untouched `interval`. Nothing downstream re-derives
    // either one from a filename.
    intervals_ch = chroms_out
        .splitText()
        .map { line -> line.trim() }
        .filter { line -> line }
        // Primary contigs only (chr1-22/X/Y/M/MT, prefix optional): alt/decoy/HLA names
        // (e.g. HLA-DQB1*02:01:01) cannot be expressed as a bcftools/pysam region string
        // -- colons parse as coordinates -- and CADD/SpliceAI/ANNOVAR carry no data for
        // them anyway.
        .filter { iv ->
            def keep = iv ==~ /(chr)?([0-9]{1,2}|X|Y|M|MT)/
            if (!keep) log.warn "annotate_snps: skipping non-primary contig '${iv}' -- not expressible as a region string and not covered by any annotation source"
            keep
        }
        .map { iv -> record(interval: iv, tag: interval_tag(iv)) }

    // vcf_tbi carries exactly one (vcf, tbi) item (main.nf wraps the two params in a
    // single channel.of(tuple(...))). Passing it as a singleton alongside the multi-item
    // intervals_ch lets Nextflow's process-call broadcast pair it against every interval
    // -- the same pattern used for squirls'/spliceai's singleton config/reference
    // arguments below.
    // Destructured before the process sees it -- see nf/shared/modules/other/tools.nf.
    iv_parts = intervals_ch.multiMap { r ->
        interval: r.interval
        tag: r.tag
    }
    split = split_vcf(vcf_tbi.first(), iv_parts.interval, iv_parts.tag)
    vcfs_ch = split.out

    // A supplied CADD table no longer skips scoring by itself: it prescores. The cadd
    // subworkflow copies its scores per interval and sends only the remainder to CADD.sh,
    // so a table that covers everything costs a lookup, and one that covers nothing is a
    // full run -- either way the merged output (prescored + newly scored) is republished to
    // 01_cadd, i.e. a valid precomputed_cadd for the next run.
    //
    // With skip_spliceai_squirls ALSO set the table passes through untouched -- the
    // ANNOVAR-only re-run, where the table is asserted complete and in the callset's contig
    // naming (main.nf warns; nothing verifies either). build_cadd_humandb still runs below
    // in every case -- reformatting and indexing that table into ANNOVAR's humandb layout is
    // a different job from scoring a genome, and cheap next to it.
    if (precomputed_cadd && skip_spliceai_squirls) {
        cadd_out = channel.value(tuple(precomputed_cadd[0], precomputed_cadd[1]))
        cadd_versions = channel.empty()
    } else {
        cadd_res = cadd(vcfs_ch, cadd_data_dir, precomputed_cadd)
        cadd_out = cadd_res.cadd_merged
        cadd_versions = cadd_res.versions
    }

    // SpliceAI and SQUIRLS skip together, never separately. The only VCF this pipeline
    // publishes is post-SQUIRLS, so there is no artifact to re-enter at the point between
    // them -- separate flags would name a resume point that does not exist.
    //
    // Skipping ASSERTS the input VCF already carries both tools' INFO fields. Nothing here
    // verifies that: a plain VCF passed with the flag set yields empty SpliceAI/SQUIRLS
    // columns rather than an error, because add_splice_scores reads absent INFO as absent.
    if (skip_spliceai_squirls) {
        annovar_vcfs = vcfs_ch
        splice_versions = channel.empty()
        // The supplied input already IS the annotated VCF; re-merging the split parts would
        // republish a copy of what the caller just handed in.
        annotated_out = channel.empty()
    } else {
        spliceai_res = spliceai(vcfs_ch, reference, spliceai_precomputed_scores, spliceai_precomputed_tbi)
        squirls_res = squirls(spliceai_res.vcf, squirls_config, squirls_jannovar)
        annovar_vcfs = squirls_res.compressed

        // The per-interval VCFs recombined into one cohort VCF, published so a later ANNOVAR
        // re-run can start from here instead of repeating CADD/SpliceAI/SQUIRLS. Named for
        // its actual contents -- see merge_vcfs.nf.
        // Not named annotated_vcf: that is the emit name below, and an emit shadows a
        // same-named local, silently breaking `.versions` on it.
        // toSortedList, not toList: collected items arrive in completion order. That varied
        // the VCF's interval order and the list's hash -- but worse, these are TWO separate
        // collections off one channel, so nothing guaranteed the tbi list matched the vcf
        // list, pairing each VCF with another interval's index. Both share the `tag` stem,
        // so sorting on the filename fixes both.
        merged_vcf = merge_vcfs(
            squirls_res.compressed.map { r -> r.vcf }
                .toSortedList { a, b -> a.name <=> b.name },
            squirls_res.compressed.map { r -> r.tbi }
                .toSortedList { a, b -> a.name <=> b.name },
            "${cohort}_${data_type}".toString()
        )
        annotated_out = merged_vcf.out
        splice_versions = spliceai_res.versions.mix(squirls_res.versions, merged_vcf.versions)
    }

    // Run-once, then broadcast to every per-interval ANNOVAR task. No .first() here: both
    // inputs are values, so the output is already a value channel and .first() on one only
    // earns `WARN: The operator 'first' is useless when applied to a value channel`.
    cadd_humandb = build_cadd_humandb(cadd_out, annovar_dir)

    // annovar_vcfs carries (interval, tag, vcf, tbi) -- the SpliceAI and SQUIRLS INFO fields
    // ANNOVAR's add_splice_scores step lifts into avinput columns are already in those VCFs,
    // whether this run produced them or the caller supplied them.
    annovar_res = annovar(
        annovar_vcfs,
        cadd_humandb.out,
        cohort,
        data_type,
        reference,
        clinvar_date,
        annovar_protocols,
        annovar_operations,
        omim_xref,
        annovar_dir,
        splice_scores_script
    )

    emit:
    annovar_tsv: Path = annovar_res.annovar_tsv
    // Published in its own right, not just consumed by build_cadd_humandb: the merged
    // genome-wide CADD table is expensive to regenerate and useful outside this run.
    cadd_merged: Tuple<Path, Path> = cadd_out
    // SpliceAI + SQUIRLS only. Together with cadd_merged it is everything a re-annotation
    // needs; on its own it is not a fully annotated VCF. Empty when skip_spliceai_squirls
    // is set -- the caller already holds this file, so nothing is republished.
    annotated_vcf: Tuple<Path, Path> = annotated_out
    versions: Channel<Path> = chroms.versions
        .mix(
            split.versions,
            cadd_versions,
            splice_versions,
            cadd_humandb.versions,
            annovar_res.versions
        )
}
