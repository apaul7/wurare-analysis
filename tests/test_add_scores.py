#!/usr/bin/env python3
"""Unit tests for nf/annotate_snps/bin/add_scores.py.

    python3 tests/test_add_scores.py

That script decides, per record, whether a variant already has a SpliceAI score in the
precomputed set (cheap, copy it) or has to be sent to SpliceAI proper (expensive). Getting it
wrong is silent in both directions: a false match attaches someone else's score, and a missed
match only costs compute, so neither shows up as an error.

The multi-allelic case is the one to watch. The upstream copy this was ported from EXITS on a multi-allelic
input record, because `vt decompose` runs upstream of it there and one arriving intact proves
a broken pipeline. This pipeline has no merge stage and decomposes later (`bcftools norm
-m-both` in normalize_vcf, part of the ANNOVAR chain), so multi-allelic records arrive on an
ordinary run and are routed to `unscored` instead. If someone ever syncs that hunk back from
upstream, the multi-allelic test below is what fails.

Invoked by path, not through PATH: Nextflow's bin/ is only on PATH inside a task.

Needs pysam. Skips cleanly without it.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import pysam
except ImportError:
    print("SKIP: pysam not installed")
    sys.exit(0)

SCRIPT = (Path(__file__).parent.parent / "nf" / "annotate_snps" / "bin" / "add_scores.py")

CONTIGS = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
"""
SPLICEAI_INFO = ('##INFO=<ID=SpliceAI,Number=.,Type=String,'
                 'Description="SpliceAI annotation">\n')
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

# Only the precomputed set declares INFO/SpliceAI. The input VCF must NOT: add_scores.py adds
# that line to the output header itself, and pysam rejects a duplicate id outright. A split
# VCF coming off the real cohort input has no SpliceAI field yet either, so this matches what
# the process actually receives.
INPUT_HEADER = CONTIGS + COLUMNS
PRECOMPUTED_HEADER = CONTIGS + SPLICEAI_INFO + COLUMNS

SCORE = "G|GENE|0.01|0.00|0.00|0.00|-7|-2|35|-30"

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def write_vcf(path, header, rows):
    """Writes a bgzipped, tabix-indexed VCF. `rows` are raw tab-joined record strings."""
    plain = path.with_suffix("")
    plain.write_text(header + "".join(r + "\n" for r in rows))
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def run(tmp, input_rows, precomputed_rows, region="chr1", allow_fail=False):
    """Runs add_scores.py; returns (scored_records, unscored_records, stderr).

    `allow_fail` returns (None, None, stderr) on a non-zero exit instead of raising, so a
    caller can assert the script did NOT exit. Without it, re-syncing upstream's multi-allelic
    `sys.exit` would blow this file up with a traceback instead of naming the regression.
    """
    inp = tmp / "in.vcf.gz"
    pre = tmp / "pre.vcf.gz"
    write_vcf(inp, INPUT_HEADER, input_rows)
    write_vcf(pre, PRECOMPUTED_HEADER, precomputed_rows)
    scored = tmp / "out.scored.vcf.gz"
    unscored = tmp / "out.unscored.vcf.gz"

    res = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(inp), "--precomputed", str(pre),
         "--region", region, "--scored", str(scored), "--unscored", str(unscored)],
        capture_output=True, text=True)
    if res.returncode != 0:
        if allow_fail:
            return None, None, res.stderr
        raise RuntimeError(f"add_scores.py exited {res.returncode}: {res.stderr[-500:]}")

    def read(p):
        return list(pysam.VariantFile(str(p)))

    return read(scored), read(unscored), res.stderr


def main():
    check("the script exists", SCRIPT.is_file(), str(SCRIPT))
    if not SCRIPT.is_file():
        return 1

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)

        # --- the ordinary split -------------------------------------------------------
        print("a variant present in the precomputed set is scored; one absent is not")
        scored, unscored, _err = run(
            tmp,
            input_rows=["chr1\t100\t.\tA\tG\t.\t.\t.",
                        "chr1\t200\t.\tC\tT\t.\t.\t."],
            precomputed_rows=[f"chr1\t100\t.\tA\tG\t.\t.\tSpliceAI={SCORE}"])
        check("one record scored", len(scored) == 1, f"{len(scored)}")
        check("one record unscored", len(unscored) == 1, f"{len(unscored)}")
        if scored:
            check("the score is carried over",
                  SCORE in str(scored[0].info["SpliceAI"][0]), str(scored[0].info.items()))
            check("and onto the right variant", scored[0].pos == 100, str(scored[0].pos))
        if unscored:
            check("the unmatched variant is the other one",
                  unscored[0].pos == 200, str(unscored[0].pos))

        # --- position matches but the allele does not ----------------------------------
        # The key is (chrom, pos, ref, alt); matching on position alone would attach a score
        # computed for a different substitution.
        print("same position, different ALT -- must not match")
        scored, unscored, _err = run(
            tmp,
            input_rows=["chr1\t100\t.\tA\tT\t.\t.\t."],
            precomputed_rows=[f"chr1\t100\t.\tA\tG\t.\t.\tSpliceAI={SCORE}"])
        check("nothing is scored", len(scored) == 0, f"{len(scored)}")
        check("the record goes to unscored", len(unscored) == 1, f"{len(unscored)}")

        # --- multi-allelic input: THE divergence from upstream --------------------
        # Upstream exits here. This pipeline must not: normalize_vcf decomposes later, so
        # these arrive on an ordinary run. They go to SpliceAI proper rather than being
        # matched on alts[0], which would attach a one-ALT score to a multi-ALT record.
        print("a multi-allelic record is routed to SpliceAI, not matched and not fatal")
        scored, unscored, err = run(
            tmp,
            input_rows=["chr1\t100\t.\tA\tG,T\t.\t.\t."],
            precomputed_rows=[f"chr1\t100\t.\tA\tG\t.\t.\tSpliceAI={SCORE}"],
            allow_fail=True)
        check("the script does not exit",
              scored is not None,
              f"it exited -- upstream's `vt decompose` assumption does not hold here: "
              f"{err[-300:]}")
        if scored is not None:
            check("nothing is scored despite the first ALT matching",
                  len(scored) == 0,
                  f"{len(scored)} -- a score for one ALT was attached to a multi-ALT record")
            check("the record is sent to unscored", len(unscored) == 1, f"{len(unscored)}")
            check("and the count is reported on stderr",
                  "multi-allelic" in err, err[-300:])

        # --- contig absent from the precomputed file ------------------------------------
        # pysam raises ValueError on a fetch for an unknown contig; the script treats it as
        # "nothing here is prescored" rather than dying.
        print("a contig missing from the precomputed set is not an error")
        scored, unscored, _err = run(
            tmp,
            input_rows=["chr2\t100\t.\tA\tG\t.\t.\t."],
            precomputed_rows=[f"chr1\t100\t.\tA\tG\t.\t.\tSpliceAI={SCORE}"],
            region="chr2")
        check("nothing is scored", len(scored) == 0, f"{len(scored)}")
        check("the record survives to unscored", len(unscored) == 1, f"{len(unscored)}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("add_scores.py routes records correctly")
    return 0


if __name__ == "__main__":
    sys.exit(main())
