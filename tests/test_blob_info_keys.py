#!/usr/bin/env python3
"""Regression tests for the SVDB blob-key lister.

Covers nf/annotate_svs/assets/list_blob_info_keys.awk directly -- no Nextflow, no container.

The output of this script is fed straight to `bcftools annotate -x`, which is a REMOVAL, so
both directions of a wrong answer are silent and bad. Naming too much strips real annotation
off the VCF AnnotSV reads, and the missing column just looks like AnnotSV not reporting it.
Naming too little leaves the report as wide and as unopenable as it was before.

The suffix set here is the same one promote_caller_support.awk uses to recognise a blob. If
one moves the other has to move with it, so the suffix test is the one to read first.

    python tests/test_blob_info_keys.py

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

AWK = (Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
       / "list_blob_info_keys.awk")

# The keys a real merged cohort VCF declares alongside the blobs. Every one of these has to
# survive, and SVLEN especially: it is what AnnotSV sizes a variant from.
KEPT = [
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of SV">',
    '##INFO=<ID=CALLER_SUPP,Number=.,Type=String,Description="Callers supporting this record">',
    '##INFO=<ID=NCALLER,Number=1,Type=Integer,Description="Number of supporting callers">',
    '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Support vector">',
    '##INFO=<ID=gnomad_sv_AF,Number=1,Type=Float,Description="gnomAD SV allele frequency">',
    '##INFO=<ID=PRPOS,Number=.,Type=String,Description="Breakpoint probability distribution">',
]

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run(lines):
    """Run the lister over `lines` and return its output as a list of INFO/KEY strings."""
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "header.txt"
        src.write_text("\n".join(lines) + "\n")
        r = subprocess.run(["awk", "-f", str(AWK), str(src)],
                           capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"awk failed: {r.stderr}")
    out = r.stdout.strip()
    return out.split(",") if out else []


def blob(tag, suffix):
    return (f'##INFO=<ID={tag}_{suffix},Number=1,Type=String,'
            f'Description="{suffix} of {tag}">')


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("a merged trio -- the blobs go, everything else stays")
    suffixes = ["CHROM", "POS", "QUAL", "FILTERS", "INFO", "SAMPLE", "FORMAT"]
    header = ["##fileformat=VCFv4.2"] + KEPT
    for tag in ("SAMP1", "SAMP2", "SAMP3"):
        header += [blob(tag, s) for s in suffixes]
    header += ["##contig=<ID=chr1,length=248956422>",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1"]
    keys = run(header)
    check("names every blob of every input", len(keys) == 21, f"got {len(keys)}: {keys}")
    check("prefixes each with INFO/ for bcftools -x",
          all(k.startswith("INFO/") for k in keys), keys)
    for line in KEPT:
        name = line.split("ID=")[1].split(",")[0]
        check(f"leaves {name} alone", f"INFO/{name}" not in keys, keys)

    print("the suffix set, one at a time -- keep in step with promote_caller_support.awk")
    for s in suffixes:
        check(f"_{s} is a blob", run([blob("T", s)]) == [f"INFO/T_{s}"])
    for s in ("AF", "OCC", "END", "LEN"):
        check(f"_{s} is not", run([blob("T", s)]) == [], f"stripping T_{s} would lose data")

    print("nothing to strip is a normal answer, not an error")
    check("a single-input cohort yields no keys", run(["##fileformat=VCFv4.2"] + KEPT) == [])
    check("an empty header yields no keys", run(["##fileformat=VCFv4.2"]) == [])

    print("the header line itself is parsed, not pattern-matched loosely")
    # A Description holding '<' or ',' would split into extra fields; only field 2 is read.
    tricky = ('##INFO=<ID=SAMP1_INFO,Number=1,Type=String,'
              'Description="INFO of SAMP1, <untouched>">')
    check("a comma or angle bracket in Description is harmless",
          run([tricky]) == ["INFO/SAMP1_INFO"])
    # FORMAT header lines declare per-sample fields and are a different namespace entirely.
    fmt = '##FORMAT=<ID=SAMP1_INFO,Number=1,Type=String,Description="not an INFO key">'
    check("a FORMAT line of the same name is not an INFO key", run([fmt]) == [])
    # The lister is handed a header, but being handed a whole VCF must not corrupt the list.
    rec = "chr1\t100\tsv_1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SAMP9_INFO=x"
    check("a record line contributes nothing", run([blob("T", "INFO"), rec]) == ["INFO/T_INFO"])

    print()
    if failures:
        print(f"{len(failures)} failed: {', '.join(failures)}")
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
