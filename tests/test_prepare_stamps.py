#!/usr/bin/env python3
"""Regression tests for the Phase 1 ID and ALGORITHMS stamps.

Covers nf/annotate_svs/assets/stamp_records.awk directly -- no Nextflow, no container, so
this runs in under a second and can gate every change to the stamping logic.

The bug that prompted these: an input already carrying INFO/ALGORITHMS kept its stale
value instead of being restamped from the samplesheet's caller, so a VCF that had passed
through another pipeline silently mislabelled the entire callset. Nothing failed; the
output just said the wrong thing.

    python tests/test_prepare_stamps.py

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

AWK = Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets" / "stamp_records.awk"

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1"""

ALG_HEADER = ('##INFO=<ID=ALGORITHMS,Number=.,Type=String,'
              'Description="Source algorithms">')

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def stamp(body_lines, tag="SET.tiddit", caller="tiddit", extra_header=None):
    """Run the awk over a synthetic VCF, return (header_lines, record_fields)."""
    header = HEADER if extra_header is None else HEADER.replace(
        "#CHROM", extra_header + "\n#CHROM")
    vcf = header + "\n" + "\n".join(body_lines) + "\n"
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.vcf"
        src.write_text(vcf)
        res = subprocess.run(
            ["awk", "-v", f"TAG={tag}", "-v", f"CALLER={caller}", "-f", str(AWK), str(src)],
            capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"awk failed: {res.stderr.strip()[:300]}")
    lines = res.stdout.splitlines()
    return ([l for l in lines if l.startswith("#")],
            [l.split("\t") for l in lines if not l.startswith("#")])


def rec(pos, info, rid="."):
    return f"chr1\t{pos}\t{rid}\tN\t<DEL>\t99\tPASS\t{info}\tGT\t0/1"


def info_of(record):
    return dict(kv.split("=", 1) for kv in record[7].split(";") if "=" in kv)


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("ALGORITHMS is stamped from the caller, never inherited")
    # THE regression: input already says manta, samplesheet says tiddit. tiddit must win.
    _, recs = stamp([rec(1000, "SVTYPE=DEL;END=2000;SVLEN=-1000;ALGORITHMS=manta")])
    check("a pre-existing ALGORITHMS is overridden",
          info_of(recs[0])["ALGORITHMS"] == "tiddit",
          f"got {info_of(recs[0]).get('ALGORITHMS')!r} -- a stale value silently "
          "mislabels the whole callset")
    check("no stale value survives anywhere in INFO",
          "manta" not in recs[0][7], f"INFO={recs[0][7]}")

    # Position in the INFO string must not matter -- callers put it anywhere.
    for label, info in [("first", "ALGORITHMS=manta;SVTYPE=DEL;END=2000"),
                        ("middle", "SVTYPE=DEL;ALGORITHMS=manta;END=2000"),
                        ("last", "SVTYPE=DEL;END=2000;ALGORITHMS=manta")]:
        _, r = stamp([rec(1000, info)])
        check(f"overridden when ALGORITHMS is {label} in INFO",
              info_of(r[0])["ALGORITHMS"] == "tiddit" and "manta" not in r[0][7],
              f"INFO={r[0][7]}")

    # Other INFO keys are collateral if the strip is written carelessly.
    _, r = stamp([rec(1000, "SVTYPE=DEL;ALGORITHMS=manta;END=2000;SVLEN=-1000")])
    got = info_of(r[0])
    check("neighbouring INFO keys survive the strip",
          got.get("SVTYPE") == "DEL" and got.get("END") == "2000"
          and got.get("SVLEN") == "-1000", f"INFO={r[0][7]}")

    _, r = stamp([rec(1000, ".")])
    check("empty INFO gets the stamps without a leading ';'",
          not r[0][7].startswith(";") and r[0][7].split(";")[0] == "ALGORITHMS=tiddit",
          f"INFO={r[0][7]!r}")

    print("CALLER_SUPP/NCALLER are stamped for records that bypass axis A")
    # Without this, the field is present on merged records and absent on unmerged ones, and
    # a filter written against it means different things on different rows.
    _, r = stamp([rec(1000, "SVTYPE=DEL")])
    got = info_of(r[0])
    check("CALLER_SUPP is stamped from the caller", got.get("CALLER_SUPP") == "tiddit",
          f"INFO={r[0][7]}")
    check("NCALLER starts at 1", got.get("NCALLER") == "1", f"INFO={r[0][7]}")

    _, r = stamp([rec(1000, "SVTYPE=DEL;CALLER_SUPP=manta,delly;NCALLER=2")])
    got = info_of(r[0])
    check("a stale CALLER_SUPP from a previous run is replaced",
          got.get("CALLER_SUPP") == "tiddit" and got.get("NCALLER") == "1",
          f"INFO={r[0][7]}")
    check("no duplicate CALLER_SUPP is left behind",
          r[0][7].count("CALLER_SUPP=") == 1, f"INFO={r[0][7]}")

    _, r = stamp([rec(1000, "SVTYPE=DEL;END=2000")])
    check("INFO with no prior ALGORITHMS gains it",
          info_of(r[0])["ALGORITHMS"] == "tiddit" and "SVTYPE=DEL" in r[0][7],
          f"INFO={r[0][7]}")

    print("record IDs are unique and stable")
    _, recs = stamp([rec(1000, "SVTYPE=DEL"), rec(1000, "SVTYPE=DUP"),
                     rec(2000, "SVTYPE=DEL")])
    ids = [r[2] for r in recs]
    check("two records at one position get distinct IDs", len(set(ids)) == 3,
          f"ids={ids} -- a CHROM_POS scheme collides here, which is why it is a counter")
    check("IDs are TAG-prefixed and sequential",
          ids == ["SET.tiddit_1", "SET.tiddit_2", "SET.tiddit_3"], f"ids={ids}")

    _, recs = stamp([rec(1000, "SVTYPE=DEL", rid="manta_1")])
    check("a pre-existing ID is replaced", recs[0][2] == "SET.tiddit_1",
          f"got {recs[0][2]!r}")

    print("header handling")
    hdr, _ = stamp([rec(1000, "SVTYPE=DEL")])
    check("ALGORITHMS header is added when absent",
          sum(1 for l in hdr if l.startswith("##INFO=<ID=ALGORITHMS")) == 1,
          "downstream bcftools needs the key declared to read it typed")

    hdr, _ = stamp([rec(1000, "SVTYPE=DEL;ALGORITHMS=manta")], extra_header=ALG_HEADER)
    check("ALGORITHMS header is not duplicated when already declared",
          sum(1 for l in hdr if l.startswith("##INFO=<ID=ALGORITHMS")) == 1,
          "a duplicate declaration is malformed VCF")

    hdr, recs = stamp([rec(1000, "SVTYPE=DEL")])
    check("the #CHROM line is preserved exactly once",
          sum(1 for l in hdr if l.startswith("#CHROM")) == 1)
    check("record count is unchanged", len(recs) == 1, f"got {len(recs)}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("all stamps behave")
    return 0


if __name__ == "__main__":
    sys.exit(main())
