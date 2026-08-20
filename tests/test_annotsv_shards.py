#!/usr/bin/env python3
"""Regression tests for the AnnotSV input sharding and output reassembly.

Covers nf/annotate_svs/assets/shard_vcf.awk and concat_annotsv_tsv.awk directly -- no
Nextflow, no container, no AnnotSV.

    python tests/test_annotsv_shards.py

WHY THESE ARE UNIT TESTS. AnnotSV itself cannot be exercised here: annotate_cohort runs only
when all four annotation params are set, and the annotation bundle is a multi-GB site-local
install, so the whole branch is skipped in every e2e suite. The sharding logic is therefore
pushed into two standalone awk assets precisely so the parts that can silently corrupt the
report are testable without it. What is left untestable here -- whether AnnotSV annotates a
record the same way in a shard as in the whole cohort -- is checked on real data by the
published annotsv_coverage.txt.

WHAT SHARDING IS FOR. AnnotSV is Tcl, and Tcl 8 caps a single value at 2^31-1 bytes. A
300-sample cohort VCF is ~9.2 GB after stripping, so AnnotSV dies in its first step with
"max size for a Tcl value (2147483647 bytes) exceeded". Splitting the input is what gets
under that without dropping the genotype columns.

The two failure modes worth guarding, both silent:

  * a record lost or duplicated by the split -- the report is quietly incomplete
  * a shard header surviving into the middle of the concatenated TSV -- both consumers read
    their columns BY NAME from line 1, so a stray header becomes a data row and every
    column read after it can come from the wrong place
"""

import subprocess
import sys
import tempfile
from pathlib import Path

ASSETS = Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
SHARD_AWK = ASSETS / "shard_vcf.awk"
CONCAT_AWK = ASSETS / "concat_annotsv_tsv.awk"

HEADER = ("##fileformat=VCFv4.2\n"
          "##contig=<ID=chr1,length=248956422>\n"
          "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1")

TSV_HEADER = "AnnotSV_ID\tID\tAnnotation_mode\tACMG_class"

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def body(n):
    """n records, ascending POS, as stamp_records.awk would mint them.

    IDs are zero-padded so every record is the same byte length. The split is byte-budgeted,
    so equal-width records make the shard boundaries exact and the assertions below readable.
    """
    return [f"chr1\t{1000 + i * 10}\tSET.manta_{i:03d}\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL\tGT\t0/1"
            for i in range(1, n + 1)]


def budget_for(records, per_shard, header=HEADER):
    """Byte budget that fits exactly `per_shard` equal-width records beside the header."""
    return len(header) + 1 + per_shard * (len(records[0]) + 1)


def shard(tmp, records, max_bytes, header=HEADER):
    """Runs shard_vcf.awk. Returns (CompletedProcess, sorted list of shard Paths)."""
    hdr = tmp / "header.vcf"
    hdr.write_text(header + "\n")
    res = subprocess.run(["awk", "-v", f"MAXBYTES={max_bytes}", "-v", f"HDR={hdr}",
                          "-f", str(SHARD_AWK)],
                         input="\n".join(records) + ("\n" if records else ""),
                         capture_output=True, text=True, cwd=str(tmp))
    return res, sorted(tmp.glob("shard_*.vcf"))


def concat(tmp, shard_texts):
    """Writes each string as a shard TSV and concatenates them. Returns CompletedProcess."""
    paths = []
    for i, text in enumerate(shard_texts, start=1):
        p = tmp / f"shard_{i:05d}.annotsv.tsv"
        p.write_text(text)
        paths.append(str(p))
    return subprocess.run(["awk", "-f", str(CONCAT_AWK)] + paths,
                          capture_output=True, text=True)


def main():
    for asset in (SHARD_AWK, CONCAT_AWK):
        if not asset.is_file():
            sys.exit(f"missing {asset}")

    print("the split loses no records and duplicates none")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        records = body(25)
        budget = budget_for(records, 10)
        res, shards = shard(tmp, records, budget)
        check("awk succeeds", res.returncode == 0, res.stderr.strip()[:200])
        check("25 records at a 10-record budget makes 3 shards", len(shards) == 3,
              f"got {[p.name for p in shards]}")
        check("and the count is reported on stdout", res.stdout.strip() == "3",
              f"got {res.stdout.strip()!r}")

        n_hdr = len(HEADER.splitlines())
        emitted = []
        for p in shards:
            lines = p.read_text().splitlines()
            check(f"{p.name} carries the whole header",
                  lines[:n_hdr] == HEADER.splitlines(),
                  "a shard without the header is not a VCF AnnotSV can read")
            emitted += lines[n_hdr:]

        check("every record survives, exactly once, in order", emitted == records,
              f"{len(emitted)} records out vs {len(records)} in")
        check("the shards are 10/10/5",
              [len(p.read_text().splitlines()) - n_hdr for p in shards] == [10, 10, 5])

        # The budget is the whole point, so it is the thing asserted -- including the header,
        # which is repeated in every shard and counts against every shard.
        sizes = [p.stat().st_size for p in shards]
        check("no shard exceeds the byte budget", all(s <= budget for s in sizes),
              f"budget={budget} sizes={sizes}")

    print("the budget is honoured whatever the record width")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # The defect this replaces: a record-count budget cannot know how wide a record is.
        # A cohort with 300 sample columns spends ~12.6 KB per record where a 1-sample cohort
        # spends ~60 bytes, and the count that was safe for one blew the ceiling on the other.
        wide = [r + "\t0/1" * 300 for r in body(9)]
        budget = budget_for(wide, 2)
        res, shards = shard(tmp, wide, budget)
        check("wide records still respect the budget",
              all(p.stat().st_size <= budget for p in shards),
              f"budget={budget} sizes={[p.stat().st_size for p in shards]}")
        check("and split into more shards, not bigger ones", len(shards) == 5,
              f"got {len(shards)}")
        n_hdr = len(HEADER.splitlines())
        out = []
        for p in shards:
            out += p.read_text().splitlines()[n_hdr:]
        check("with every wide record preserved", out == wide, f"{len(out)} of {len(wide)}")

        # A record that cannot fit a shard of its own is unsplittable; there is no smaller
        # unit than a VCF record. Better to say so than to emit an over-budget shard.
        res, _ = shard(tmp, body(3), len(HEADER) + 1 + 5)
        check("a record too big for any shard is fatal", res.returncode != 0,
              f"rc={res.returncode}")
        check("and says which record and what would fit",
              "does not fit" in res.stderr, res.stderr.strip()[:160])

    print("the split refuses inputs it cannot shard correctly")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # An unreadable header would yield shards that are not VCFs. AnnotSV would reject
        # them one at a time without ever saying the header was the problem.
        res = subprocess.run(["awk", "-v", "MAXBYTES=100000", "-v", f"HDR={tmp}/nope.vcf",
                              "-f", str(SHARD_AWK)],
                             input="\n".join(body(3)) + "\n",
                             capture_output=True, text=True, cwd=str(tmp))
        check("a missing header file is fatal", res.returncode != 0, f"rc={res.returncode}")
        check("and says so", "header" in res.stderr.lower(), res.stderr.strip()[:120])

        res, _ = shard(tmp, [], 100000)
        check("an empty body is fatal rather than silently producing nothing",
              res.returncode != 0, f"rc={res.returncode}")

    print("the concatenation keeps exactly one header, first")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        res = concat(tmp, [
            f"{TSV_HEADER}\n1_1000_2000_DEL_1\tSET.manta_1\tsplit\t4\n",
            f"{TSV_HEADER}\n1_3000_4000_DEL_1\tSET.manta_2\tsplit\t3\n",
        ])
        check("awk succeeds", res.returncode == 0, res.stderr.strip()[:200])
        out = res.stdout.splitlines()
        check("the header is line 1", out[0] == TSV_HEADER, f"got {out[0]!r}")
        check("and appears exactly once", out.count(TSV_HEADER) == 1,
              "a repeated header becomes a data row whose cells are column names")
        check("every shard's rows are present, in shard order",
              [l.split("\t")[1] for l in out[1:]] == ["SET.manta_1", "SET.manta_2"],
              f"got {out[1:]}")

    print("degenerate shards are tolerated, inconsistent ones are not")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        res = concat(tmp, [
            f"{TSV_HEADER}\n1_1000_2000_DEL_1\tSET.manta_1\tsplit\t4\n",
            "",                       # AnnotSV produced nothing for this shard
            f"{TSV_HEADER}\n",        # header only: no record here could be typed
            f"{TSV_HEADER}\n1_5000_6000_DEL_1\tSET.manta_3\tsplit\t5\n",
        ])
        out = res.stdout.splitlines()
        check("empty and header-only shards are not an error", res.returncode == 0,
              res.stderr.strip()[:200])
        check("they contribute no rows", out.count(TSV_HEADER) == 1 and len(out) == 3,
              f"got {out}")

        # AnnotSV's column set shifts between releases and bundles. Two shards annotated
        # differently, concatenated under one header, misalign every later column silently.
        res = concat(tmp, [
            f"{TSV_HEADER}\n1_1000_2000_DEL_1\tSET.manta_1\tsplit\t4\n",
            f"{TSV_HEADER}\tExtra_Column\n1_3000_4000_DEL_1\tSET.manta_2\tsplit\t3\tx\n",
        ])
        check("a shard with a different header is fatal", res.returncode != 0,
              f"rc={res.returncode}")
        check("and names the offending shard", "shard_00002" in res.stderr,
              res.stderr.strip()[:160])

        res = concat(tmp, ["", ""])
        check("all-empty input is fatal rather than an empty report",
              res.returncode != 0, f"rc={res.returncode}")

    print("split then reassemble is lossless end to end")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        records = body(23)
        res, shards = shard(tmp, records, budget_for(records, 5))
        check("5 shards for 23 records", len(shards) == 5, f"got {len(shards)}")

        # Stand in for AnnotSV: one TSV row per input record, carrying its ID.
        fake = []
        for p in shards:
            rows = [l for l in p.read_text().splitlines() if not l.startswith("#")]
            fake.append(TSV_HEADER + "\n"
                        + "".join(f"x\t{r.split(chr(9))[2]}\tsplit\t3\n" for r in rows))
        res = concat(tmp, fake)
        ids = [l.split("\t")[1] for l in res.stdout.splitlines()[1:]]
        check("every record reaches the assembled TSV exactly once",
              ids == [f"SET.manta_{i:03d}" for i in range(1, 24)],
              f"{len(ids)} ids out of 23")

    print()
    if failures:
        print(f"{len(failures)} FAILED: " + ", ".join(failures))
        return 1
    print("sharding and reassembly hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
