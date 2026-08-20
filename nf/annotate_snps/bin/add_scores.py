#!/usr/bin/env python3
"""Split a VCF into records that already have a precomputed SpliceAI score and records that do not.

Called by the `add_precomputed` process in modules/annotations/spliceai.nf. Everything written
to `--scored` carries an INFO/SpliceAI field copied from the precomputed set; everything
written to `--unscored` goes on to run SpliceAI proper, which is the expensive step this
lookup exists to avoid.

Lives in `bin/` rather than inline in the module: Nextflow puts `$projectDir/bin` on PATH
automatically, so the process can call `add_scores.py` directly. As a heredoc this had no
syntax checking, no linting and no way to test it.

Note `bin/` is resolved against the *entry script's* project directory, so this is on PATH for
`nextflow run nf/annotate_snps/main.nf` but not for a harness launched from elsewhere. A test
should invoke this file by path instead of relying on PATH.

Adapted from an upstream pipeline. ONE DELIBERATE DIVERGENCE, in how multi-allelic input
records are handled -- see the comment in the record loop. Do not sync that hunk back and
forth without reading it; the two pipelines decompose at different points.
"""

import argparse
import sys

import pysam


def variant_key(record):
    """(chrom, pos, ref, alt) for a biallelic record."""
    return (record.chrom, record.pos, record.ref, record.alts[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--precomputed", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--scored", required=True)
    parser.add_argument("--unscored", required=True)
    args = parser.parse_args()

    print(f"in: {args.input}")
    print(f"region: {args.region}")
    print(f"scored out: {args.scored}")
    print(f"unscored out: {args.unscored}")

    # Region-restricted. The precomputed scores file is genome-wide, but this task only
    # ever looks up variants from one interval, so scanning the whole file in every
    # per-interval task was O(intervals x genome) for no additional matches: a lookup
    # outside this region can never be hit by the input VCF below.
    precomputed_vcf = pysam.VariantFile(args.precomputed)
    try:
        region_records = precomputed_vcf.fetch(region=args.region)
    except ValueError:
        # Contig absent from the precomputed file: nothing here is prescored. The
        # full-file scan this replaced reached the same answer by simply finding no match.
        region_records = []

    scored_lookup = {}
    n_precomputed_multiallelic = 0
    for record in region_records:
        # The lookup key uses alts[0], so any further ALTs on a precomputed record are not
        # reachable. Counted and reported rather than passed over in silence -- the matching
        # behaviour is unchanged, but a precomputed set that would benefit from being
        # decomposed no longer looks identical to one that needs nothing.
        if record.alts is not None and len(record.alts) > 1:
            n_precomputed_multiallelic += 1
        scored_lookup[variant_key(record)] = record.info["SpliceAI"]

    if n_precomputed_multiallelic:
        print(
            f"warning: {n_precomputed_multiallelic} multi-allelic record(s) in the precomputed "
            f"set; only their first ALT can be matched",
            file=sys.stderr,
        )

    input_vcf = pysam.VariantFile(args.input)

    header = input_vcf.header.copy()
    header.info.add(
        "SpliceAI",
        number=".",
        type="String",
        description="SpliceAIv1.3.1 variant annotation. These include delta scores and delta positions for acceptor gain, acceptor loss, donor gain, and donor loss."
    )

    unscored_out = pysam.VariantFile(args.unscored, "w", header=input_vcf.header)
    scored_out = pysam.VariantFile(args.scored, "w", header=header)

    n_scored = 0
    n_unscored = 0
    n_multiallelic = 0

    for record in input_vcf:
        # DIVERGENCE FROM the upstream copy, which exits with an error here instead.
        #
        # There, every record has been through `vt decompose` in merge_variants before it
        # arrives, so a multi-allelic record proves an upstream step did not run and exiting
        # is right. This pipeline has no merge stage: --vcf is whatever cohort VCF the caller
        # handed in, and the only decomposition (`bcftools norm -m-both` in normalize_vcf) is
        # part of the ANNOVAR chain, which runs AFTER this. Multi-allelic records therefore
        # reach this script on a perfectly ordinary run, and exiting would abort it.
        #
        # They are routed to `unscored` without consulting the lookup. variant_key sees only
        # alts[0], so a match would attach a score derived from one ALT to a record carrying
        # several -- a wrong annotation rather than a missing one. Sending them to SpliceAI
        # proper costs compute and returns a score per ALT, which is the correct answer.
        n_alts = 0 if record.alts is None else len(record.alts)
        if n_alts > 1:
            n_multiallelic += 1
            unscored_out.write(record)
            n_unscored += 1
            continue

        key = variant_key(record)
        if key in scored_lookup:
            new_record = record.copy()
            new_record.translate(header)
            new_record.info["SpliceAI"] = scored_lookup[key]
            scored_out.write(new_record)
            n_scored += 1
        else:
            unscored_out.write(record)
            n_unscored += 1

    if n_multiallelic:
        print(
            f"warning: {n_multiallelic} multi-allelic record(s) in {args.input} were sent to "
            f"SpliceAI rather than matched against the precomputed set",
            file=sys.stderr,
        )

    print(f"scored: {n_scored} variants")
    print(f"unscored: {n_unscored} variants")


if __name__ == "__main__":
    main()
