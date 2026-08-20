#!/usr/bin/env python3
"""Split a VCF into records already scored in a prescored CADD table and records that are not.

Called by the `extract_prescored` process in modules/annotations/cadd.nf. Every match is
written to `--scored-tsv` as a CADD output row copied from the prescored table; everything
else goes to `--unscored` and on to CADD.sh proper, which is the expensive step this lookup
exists to avoid.

The prescored table is `#Chrom Pos Ref Alt RawScore PHRED`, bgzipped and tabix-indexed
(`tabix -s 1 -b 2 -e 2`). Its contigs may use either naming this pipeline sees: the
callset's (a previous run's 01_cadd/ output, e.g. chr1/chrM) or CADD's own Ensembl naming
(raw CADD.sh output, e.g. 1/MT). The region lookup tries both. Output rows are always
written in CADD-native naming, whatever the table used, because merge_cadd re-adds the
callset's prefix to every fragment uniformly -- a fragment that kept a chr prefix would come
out as chrchr1.

Lives in `bin/` rather than inline in the module: Nextflow puts `$projectDir/bin` on PATH
automatically, so the process can call `add_cadd_scores.py` directly. Note `bin/` is
resolved against the *entry script's* project directory, so a test should invoke this file
by path instead of relying on PATH.

Structured after add_scores.py (the SpliceAI equivalent), including its multi-allelic
stance: this pipeline decomposes late (`bcftools norm -m-both` inside the ANNOVAR chain),
so multi-allelic records arrive here on an ordinary run and are routed to `unscored` --
a first-ALT match would attach one ALT's score to a record carrying several, and CADD
scores every ALT itself.
"""

import argparse
import sys

import pysam


def native_contig(contig):
    """CADD v1.6 GRCh38 naming: no chr prefix, and the mitochondrion is MT, not M."""
    c = contig[3:] if contig.startswith("chr") else contig
    return "MT" if c == "M" else c


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--prescored", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--scored-tsv", required=True)
    parser.add_argument("--unscored", required=True)
    args = parser.parse_args()

    print(f"in: {args.input}")
    print(f"region: {args.region}")
    print(f"scored out: {args.scored_tsv}")
    print(f"unscored out: {args.unscored}")

    # Region-restricted, like add_scores.py -- a lookup outside this interval can never be
    # hit by the input VCF below. The table's own contig naming is unknown, so the region is
    # resolved against what the index actually contains: as given, chr-stripped, then M->MT.
    tbx = pysam.TabixFile(args.prescored)
    # Sliced, not str.removeprefix: the analysis container's python predates 3.9.
    stripped = args.region[3:] if args.region.startswith("chr") else args.region
    candidates = []
    for c in (args.region, stripped, native_contig(args.region)):
        if c not in candidates:
            candidates.append(c)
    tsv_contig = next((c for c in candidates if c in tbx.contigs), None)

    scored_lookup = {}
    if tsv_contig is None:
        # Contig absent from the table under every naming: nothing here is prescored. Same
        # stance as add_scores.py's fetch ValueError -- an empty lookup, not an error.
        print(f"none of {candidates} in the prescored table; nothing here is prescored")
    else:
        for line in tbx.fetch(tsv_contig):
            chrom, pos, ref, alt, raw, phred = line.rstrip("\n").split("\t")[:6]
            scored_lookup[(native_contig(chrom), int(pos), ref, alt)] = (raw, phred)

    # The fragment's header comes from the table itself, so the version merge_cadd's first
    # file carries (and build_cadd_humandb parses) is the one the scores were made with.
    header_lines = list(tbx.header)
    if not any(l.startswith("##CADD") for l in header_lines):
        header_lines = ["##CADD GRCh38-v1.6 (header synthesized by add_cadd_scores.py)",
                        "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED"]

    input_vcf = pysam.VariantFile(args.input)
    unscored_out = pysam.VariantFile(args.unscored, "wz", header=input_vcf.header)

    n_scored = 0
    n_unscored = 0
    n_multiallelic = 0

    with open(args.scored_tsv, "w") as scored_out:
        for l in header_lines:
            scored_out.write(l + "\n")

        for record in input_vcf:
            # Routed to CADD proper rather than matched on alts[0] -- see the module
            # docstring, and add_scores.py's identically-reasoned hunk.
            n_alts = 0 if record.alts is None else len(record.alts)
            if n_alts > 1:
                n_multiallelic += 1
                unscored_out.write(record)
                n_unscored += 1
                continue

            key = (native_contig(record.chrom), record.pos, record.ref, record.alts[0])
            if key in scored_lookup:
                raw, phred = scored_lookup[key]
                scored_out.write("\t".join((key[0], str(record.pos), record.ref,
                                            record.alts[0], raw, phred)) + "\n")
                n_scored += 1
            else:
                unscored_out.write(record)
                n_unscored += 1

    unscored_out.close()

    if n_multiallelic:
        print(
            f"warning: {n_multiallelic} multi-allelic record(s) in {args.input} were sent to "
            f"CADD rather than matched against the prescored table",
            file=sys.stderr,
        )

    print(f"scored: {n_scored} variants")
    print(f"unscored: {n_unscored} variants")


if __name__ == "__main__":
    main()
