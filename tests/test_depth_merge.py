#!/usr/bin/env python3
"""Regression tests for the duphold FORMAT recombination.

Covers nf/annotate_svs/assets/merge_depth.awk directly -- no Nextflow, no container.

duphold takes one alignment at a time, so the cohort VCF is split per sample and put back
together. Every failure mode here is silent and severe: depth
attached to the wrong sample or the wrong record does not error, it just quietly drives the
depth filter from the wrong numbers. The join is on the unique record ID stamped in Phase 1,
never on position and never on POS/REF/ALT -- two symbolic-ALT records can share
CHROM/POS/REF/ALT and differ only in INFO/END.

    python tests/test_depth_merge.py

Needs `awk` only.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ASSETS = Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
AWK = ASSETS / "merge_depth.awk"
EXTRACT_AWK = ASSETS / "extract_depth.awk"

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1\tSAMP2\tSAMP3"""

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run_merge(records, depth_rows, header=HEADER):
    """Raw result, for the cases that are supposed to fail. Returns CompletedProcess."""
    body = "\n".join(
        f"chr1\t{pos}\t{rid}\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL\tGT\t" + "\t".join(gts)
        for rid, pos, gts in records)
    with tempfile.TemporaryDirectory() as td:
        tsv = Path(td) / "depth.tsv"
        tsv.write_text("".join("\t".join(map(str, r)) + "\n" for r in depth_rows))
        vcf = Path(td) / "cohort.vcf"
        vcf.write_text(header + "\n" + body + "\n")
        return subprocess.run(["awk", "-f", str(AWK), str(tsv), str(vcf)],
                              capture_output=True, text=True,
                              env=dict(os.environ, LC_ALL="C"))


def merge(records, depth_rows, header=HEADER):
    """records: list of (id, pos, gts). depth_rows: list of (id, sample, *values)."""
    body = "\n".join(
        f"chr1\t{pos}\t{rid}\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL\tGT\t" + "\t".join(gts)
        for rid, pos, gts in records)
    with tempfile.TemporaryDirectory() as td:
        tsv = Path(td) / "depth.tsv"
        tsv.write_text("".join("\t".join(map(str, r)) + "\n" for r in depth_rows))
        vcf = Path(td) / "cohort.vcf"
        vcf.write_text(header + "\n" + body + "\n")
        res = subprocess.run(["awk", "-f", str(AWK), str(tsv), str(vcf)],
                             capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"awk failed: {res.stderr.strip()[:300]}")
    lines = res.stdout.splitlines()
    return ([l for l in lines if l.startswith("#")],
            [l.split("\t") for l in lines if not l.startswith("#")])


def fmt_field(record, sample_idx, key, fmt_col=8):
    keys = record[fmt_col].split(":")
    vals = record[9 + sample_idx].split(":")
    return vals[keys.index(key)] if key in keys else None


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("depth lands on the right sample and the right record")
    recs = [("rec_1", 1000, ["0/1", "1/1", "0/0"]),
            ("rec_2", 2000, ["0/0", "0/1", "1/1"])]
    depth = [("rec_1", "SAMP1", 0.1, 0.2, 0.3, -4.0),
             ("rec_1", "SAMP2", 0.9, 0.8, 0.7, -0.1),
             ("rec_2", "SAMP3", 1.5, 1.6, 1.7, 3.2)]
    _, out = merge(recs, depth)
    check("SAMP1 gets its own record-1 value",
          fmt_field(out[0], 0, "DHFFC") == "0.1", f"got {fmt_field(out[0], 0, 'DHFFC')}")
    check("SAMP2 gets its own record-1 value",
          fmt_field(out[0], 1, "DHFFC") == "0.9", f"got {fmt_field(out[0], 1, 'DHFFC')}")
    check("SAMP3 gets its own record-2 value",
          fmt_field(out[1], 2, "DHFFC") == "1.5", f"got {fmt_field(out[1], 2, 'DHFFC')}")
    check("all four duphold keys are attached",
          [fmt_field(out[0], 0, k) for k in ("DHFFC", "DHBFC", "DHFC", "DHBZ")]
          == ["0.1", "0.2", "0.3", "-4.0"], f"FORMAT={out[0][8]} SAMP1={out[0][9]}")
    check("GT stays first in FORMAT and keeps its value",
          out[0][8].split(":")[0] == "GT" and out[0][9].split(":")[0] == "0/1",
          f"FORMAT={out[0][8]} SAMP1={out[0][9]}")

    print("samples and records without depth are filled, never dropped")
    check("a sample with no alignment gets '.' not a missing column",
          fmt_field(out[0], 2, "DHFFC") == "."
          and len(out[0][11].split(":")) == len(out[0][8].split(":")),
          f"SAMP3={out[0][11]} -- a short column silently shifts every later field")
    check("a record with no depth at all is still emitted with padding",
          fmt_field(out[1], 0, "DHFFC") == ".", f"got {out[1][9]}")
    check("every sample column has the same field count as FORMAT",
          all(len(out[r][9 + s].split(":")) == len(out[r][8].split(":"))
              for r in range(2) for s in range(3)))

    print("the join is on ID, not on position")
    # Same CHROM/POS/REF/ALT on both records: a positional or POS/REF/ALT join mis-pairs
    # them, an ID join cannot. This is the mis-attached-depth failure mode, made concrete.
    recs = [("rec_a", 1000, ["0/1", "0/0", "0/0"]),
            ("rec_b", 1000, ["0/0", "0/1", "0/0"])]
    depth = [("rec_a", "SAMP1", 0.05, 0.05, 0.05, -9.0),
             ("rec_b", "SAMP2", 2.50, 2.50, 2.50, 9.0)]
    _, out = merge(recs, depth)
    check("two records at one position keep their own values",
          fmt_field(out[0], 0, "DHFFC") == "0.05"
          and fmt_field(out[1], 1, "DHFFC") == "2.5",
          f"rec_a SAMP1={fmt_field(out[0], 0, 'DHFFC')} "
          f"rec_b SAMP2={fmt_field(out[1], 1, 'DHFFC')}")
    check("the second record did not inherit the first's depth",
          fmt_field(out[1], 0, "DHFFC") == ".", f"got {fmt_field(out[1], 0, 'DHFFC')}")

    print("sample order comes from the VCF, not from the depth table")
    # duphold runs finish in whatever order the executor schedules them.
    recs = [("rec_1", 1000, ["0/1", "0/1", "0/1"])]
    depth = [("rec_1", "SAMP3", 3.0, 3.0, 3.0, 3.0),
             ("rec_1", "SAMP1", 1.0, 1.0, 1.0, 1.0),
             ("rec_1", "SAMP2", 2.0, 2.0, 2.0, 2.0)]
    _, out = merge(recs, depth)
    check("reversed table order still maps by sample name",
          [fmt_field(out[0], i, "DHFFC") for i in range(3)] == ["1.0", "2.0", "3.0"],
          f"got {[fmt_field(out[0], i, 'DHFFC') for i in range(3)]}")

    print("structure")
    hdr, out = merge([("rec_1", 1000, ["0/1", "0/1", "0/1"])],
                     [("rec_1", "SAMP1", 1, 1, 1, 1)])
    for key in ("DHFFC", "DHBFC", "DHFC", "DHBZ"):
        check(f"{key} is declared in the header",
              sum(1 for l in hdr if f"##FORMAT=<ID={key}," in l) == 1,
              "bcftools needs it declared to read it typed")
    check("the #CHROM line is preserved with all samples",
          [l for l in hdr if l.startswith("#CHROM")][0].split("\t")[9:]
          == ["SAMP1", "SAMP2", "SAMP3"])
    check("record count is unchanged", len(out) == 1, f"got {len(out)}")

    _, out = merge([("rec_1", 1000, ["0/1", "0/1", "0/1"])], [])
    check("an empty depth table pads every sample rather than failing",
          all(fmt_field(out[0], i, "DHFFC") == "." for i in range(3)),
          "a cohort with no alignments at all must still produce a valid VCF")

    print("the merge join refuses unsorted input instead of under-joining")
    # This is a streaming merge join over two ID-sorted streams, not a lookup table -- the
    # table version needed an estimated 150 GB of RAM at 300 samples. The cost of streaming
    # is the precondition, and an unsorted input does not error on its own: it silently
    # drops depth for whatever it skipped past, and the depth filters then run on
    # padding. So the ordering is asserted in the awk, and these two cases prove it.
    r = run_merge([("rec_1", 1000, ["0/1", "0/1", "0/1"]),
                   ("rec_2", 2000, ["0/1", "0/1", "0/1"])],
                  [("rec_2", "SAMP1", 2, 2, 2, 2),
                   ("rec_1", "SAMP1", 1, 1, 1, 1)])
    check("a descending depth table is fatal", r.returncode != 0, f"rc={r.returncode}")
    check("and says which key broke the order", "not sorted by ID" in r.stderr,
          f"stderr={r.stderr.strip()[:120]}")

    # rec_2 before rec_1 in the BODY. recombine_depth sorts the body by ID for exactly this
    # reason; a coordinate-ordered body reaching the awk unsorted must not pass silently.
    r = run_merge([("rec_2", 1000, ["0/1", "0/1", "0/1"]),
                   ("rec_1", 2000, ["0/1", "0/1", "0/1"])],
                  [("rec_1", "SAMP1", 1, 1, 1, 1),
                   ("rec_2", "SAMP1", 2, 2, 2, 2)])
    check("a descending VCF body is fatal", r.returncode != 0, f"rc={r.returncode}")
    check("and names the body as the unsorted side", "body is not sorted" in r.stderr,
          f"stderr={r.stderr.strip()[:120]}")

    # The stream has to tolerate depth for records the cohort VCF does not contain -- a
    # sample's duphold output is cut from an earlier VCF, so a record dropped in between
    # leaves an orphan row. It must be skipped without swallowing the next record's rows.
    _, out = merge([("rec_2", 2000, ["0/1", "0/1", "0/1"])],
                   [("rec_1", "SAMP1", 9, 9, 9, 9),
                    ("rec_2", "SAMP1", 2, 2, 2, 2)])
    check("an orphan depth row is skipped, not misapplied",
          fmt_field(out[0], 0, "DHFFC") == "2", f"got {fmt_field(out[0], 0, 'DHFFC')}")

    print("extraction tolerates duphold tags that are not there")
    # bcftools query hard-errors on an undeclared tag ("no such tag defined in the VCF
    # header: FORMAT/DHBZ"), and duphold does not always emit all four -- DHBZ needs enough
    # GC-matched bins to exist. Reading FORMAT by name is tolerant by construction.
    def extract(fmt, vals, sample="SAMP1"):
        vcf = (f'##fileformat=VCFv4.2\n'
               f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n'
               f'chr1\t1000\trec_1\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL\t{fmt}\t{vals}\n')
        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "dh.vcf"
            src.write_text(vcf)
            r = subprocess.run(["awk", "-f", str(EXTRACT_AWK), str(src)],
                               capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"awk failed: {r.stderr.strip()[:200]}")
        return r.stdout.strip().split("\t")

    row = extract("GT:DHFFC:DHBFC:DHFC:DHBZ", "0/1:0.1:0.2:0.3:-4.0")
    check("all four tags present are extracted in order",
          row == ["rec_1", "SAMP1", "0.1", "0.2", "0.3", "-4.0"], f"got {row}")

    row = extract("GT:DHFFC:DHBFC:DHFC", "0/1:0.1:0.2:0.3")
    check("a missing DHBZ becomes '.' instead of an error",
          row == ["rec_1", "SAMP1", "0.1", "0.2", "0.3", "."], f"got {row}")

    row = extract("GT", "0/1")
    check("no duphold tags at all still yields a padded row",
          row == ["rec_1", "SAMP1", ".", ".", ".", "."], f"got {row}")

    row = extract("DHFC:GT:DHFFC", "0.3:0/1:0.1")
    check("tags are read by name, not by position",
          row == ["rec_1", "SAMP1", "0.1", ".", "0.3", "."], f"got {row}")

    row = extract("GT:DHFFC:DHBFC:DHFC:DHBZ", "0/1:0.1:0.2:0.3:-4.0", sample="OTHER")
    check("the sample name comes from the #CHROM line",
          row[1] == "OTHER", f"got {row}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("depth recombination holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
