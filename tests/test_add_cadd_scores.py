#!/usr/bin/env python3
"""Unit tests for nf/annotate_snps/bin/add_cadd_scores.py.

    python3 tests/test_add_cadd_scores.py

That script decides, per record, whether a variant already has a row in a prescored CADD
table (cheap, copy it) or has to be sent to CADD.sh proper (expensive). Getting it wrong is
silent in both directions: a false match attaches someone else's score, and a missed match
only costs compute, so neither shows up as an error.

The contig naming cases are the ones to watch. The table may be a previous run's 01_cadd/
output (callset-named: chr1, chrM) or raw CADD.sh output (Ensembl-named: 1, MT), while the
fragment it writes must ALWAYS be native-named -- merge_cadd re-adds the callset prefix to
every fragment, so a fragment that kept its chr prefix would come out as chrchr1, which
matches no contig anywhere and fails as a blank ANNOVAR column, not an error.

Invoked by path, not through PATH: Nextflow's bin/ is only on PATH inside a task.

Needs pysam. Skips cleanly without it.
"""

import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import pysam
except ImportError:
    print("SKIP: pysam not installed")
    sys.exit(0)

SCRIPT = (Path(__file__).parent.parent / "nf" / "annotate_snps" / "bin"
          / "add_cadd_scores.py")

VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chrM,length=16569>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

CADD_HEADER = ("##CADD GRCh38-v1.6 (c) University of Washington, Hudson-Alpha Institute "
               "for Biotechnology and Berlin Institute of Health 2013-2020. "
               "All rights reserved.\n"
               "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n")

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def write_tsv(path, rows):
    """Writes a bgzipped, tabix-indexed CADD table. `rows` are raw tab-joined strings."""
    plain = path.with_suffix("")
    plain.write_text(CADD_HEADER + "".join(r + "\n" for r in rows))
    pysam.tabix_compress(str(plain), str(path), force=True)
    # tabix -s 1 -b 2 -e 2, in pysam's 0-based column indices.
    pysam.tabix_index(str(path), seq_col=0, start_col=1, end_col=1, force=True)


def run(tmp, input_rows, table_rows, region="chr1", allow_fail=False):
    """Runs add_cadd_scores.py; returns (scored_tsv_body_rows, unscored_records, stderr).

    Body rows only -- the header's presence is asserted separately. `allow_fail` returns
    (None, None, stderr) on a non-zero exit instead of raising.
    """
    inp = tmp / "in.vcf"
    inp.write_text(VCF_HEADER + "".join(r + "\n" for r in input_rows))
    pre = tmp / "pre.tsv.gz"
    write_tsv(pre, table_rows)
    scored = tmp / "out.prescored.tsv"
    unscored = tmp / "out.unscored.vcf.gz"

    res = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(inp), "--prescored", str(pre),
         "--region", region, "--scored-tsv", str(scored), "--unscored", str(unscored)],
        capture_output=True, text=True)
    if res.returncode != 0:
        if allow_fail:
            return None, None, res.stderr
        raise RuntimeError(f"add_cadd_scores.py exited {res.returncode}: {res.stderr[-500:]}")

    body = [l for l in scored.read_text().splitlines() if not l.startswith("#")]
    return body, list(pysam.VariantFile(str(unscored))), res.stderr


def main():
    check("the script exists", SCRIPT.is_file(), str(SCRIPT))
    if not SCRIPT.is_file():
        return 1

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)

        # --- the ordinary split -------------------------------------------------------
        print("a variant present in the table is copied out; one absent is not")
        scored, unscored, _err = run(
            tmp,
            input_rows=["chr1\t100\t.\tA\tG\t.\t.\t.",
                        "chr1\t200\t.\tC\tT\t.\t.\t."],
            table_rows=["chr1\t100\tA\tG\t0.234016\t3.517"])
        check("one row scored", len(scored) == 1, f"{scored}")
        check("one record unscored", len(unscored) == 1, f"{len(unscored)}")
        if scored:
            check("the scores are carried over",
                  scored[0].endswith("0.234016\t3.517"), scored[0])
            check("the fragment row is native-named despite the chr-named table",
                  scored[0].split("\t")[0] == "1", scored[0])
        if unscored:
            check("the unmatched variant is the other one",
                  unscored[0].pos == 200, str(unscored[0].pos))

        # --- the table in CADD's own naming ---------------------------------------------
        # Raw CADD.sh output has Ensembl contigs (1, MT). The region lookup must find them
        # from a chr-named callset's interval.
        print("a native-named table matches a chr-named callset")
        scored, unscored, _err = run(
            tmp,
            input_rows=["chr1\t100\t.\tA\tG\t.\t.\t."],
            table_rows=["1\t100\tA\tG\t0.234016\t3.517"])
        check("the row is found", len(scored) == 1, f"{scored}")
        check("nothing is unscored", len(unscored) == 0, f"{len(unscored)}")

        print("the mitochondrion maps chrM -> MT")
        scored, unscored, _err = run(
            tmp,
            input_rows=["chrM\t50\t.\tA\tG\t.\t.\t."],
            table_rows=["MT\t50\tA\tG\t0.1\t1.5"],
            region="chrM")
        check("the row is found", len(scored) == 1, f"{scored}")
        if scored:
            check("and written as MT, not M or chrM",
                  scored[0].split("\t")[0] == "MT", scored[0])

        # --- position matches but the allele does not ----------------------------------
        print("same position, different ALT -- must not match")
        scored, unscored, _err = run(
            tmp,
            input_rows=["chr1\t100\t.\tA\tT\t.\t.\t."],
            table_rows=["chr1\t100\tA\tG\t0.234016\t3.517"])
        check("nothing is scored", len(scored) == 0, f"{scored}")
        check("the record goes to unscored", len(unscored) == 1, f"{len(unscored)}")

        # --- multi-allelic input --------------------------------------------------------
        # Same stance as add_scores.py: matched on alts[0] it would carry one ALT's score
        # onto a record with several, so it goes to CADD proper, which scores every ALT.
        print("a multi-allelic record is routed to CADD, not matched and not fatal")
        scored, unscored, err = run(
            tmp,
            input_rows=["chr1\t100\t.\tA\tG,T\t.\t.\t."],
            table_rows=["chr1\t100\tA\tG\t0.234016\t3.517"],
            allow_fail=True)
        check("the script does not exit", scored is not None, f"it exited: {err[-300:]}")
        if scored is not None:
            check("nothing is scored despite the first ALT matching",
                  len(scored) == 0, f"{scored}")
            check("the record is sent to unscored", len(unscored) == 1, f"{len(unscored)}")
            check("and the count is reported on stderr", "multi-allelic" in err, err[-300:])

        # --- contig absent from the table under every naming ----------------------------
        print("a contig missing from the table is not an error")
        scored, unscored, _err = run(
            tmp,
            input_rows=["chr2\t100\t.\tA\tG\t.\t.\t."],
            table_rows=["chr1\t100\tA\tG\t0.234016\t3.517"],
            region="chr2")
        check("nothing is scored", len(scored) == 0, f"{scored}")
        check("the record survives to unscored", len(unscored) == 1, f"{len(unscored)}")

        # --- the fragment header --------------------------------------------------------
        # merge_cadd reads the merged header from its first file and build_cadd_humandb
        # parses the CADD version out of the ##CADD line, so the fragment must carry both
        # lines -- copied from the table, whose version is the one the scores came from.
        print("the fragment carries the table's own header")
        run(tmp, input_rows=["chr1\t100\t.\tA\tG\t.\t.\t."],
            table_rows=["chr1\t100\tA\tG\t0.234016\t3.517"])
        head = (tmp / "out.prescored.tsv").read_text().splitlines()[:2]
        check("##CADD line first, column line second",
              head[0].startswith("##CADD GRCh38-v1.6") and head[1].startswith("#Chrom"),
              str(head))

        # --- the unscored VCF is bgzf ---------------------------------------------------
        # The process tabixes it, and tabix refuses plain gzip.
        print("the unscored output is bgzf-compressed")
        with gzip.open(tmp / "out.unscored.vcf.gz", "rb") as fh:
            fh.read(1)  # readable as gzip at all
        check("pysam can index it", _tabix_ok(tmp / "out.unscored.vcf.gz"))

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("add_cadd_scores.py routes records correctly")
    return 0


def _tabix_ok(path):
    try:
        pysam.tabix_index(str(path), preset="vcf", force=True)
        return True
    except OSError:
        return False


if __name__ == "__main__":
    sys.exit(main())
