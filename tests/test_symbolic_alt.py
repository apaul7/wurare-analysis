#!/usr/bin/env python3
"""Regression tests for assets/symbolic_alt.awk.

Covers the rewrite that keeps `gatk SVAnnotate` alive on real caller output. SVAnnotate
throws on a sequence-resolved ALT:

    java.lang.IllegalArgumentException: Unexpected ALT allele: TTTTTTCTTTCTTT...
    Expected breakpoint or symbolic ALT allele representing a structural variant record.

and it throws rather than skipping, so one such record ends the whole Talos tail. Manta emits
them as a matter of course, which is why this survived every symbolic-only fixture here and
failed on the first real cohort.

    python tests/test_symbolic_alt.py

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

AWK = Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets" / "symbolic_alt.awk"

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1"""

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run(records):
    body = "".join(f"chr1\t{p}\t{i}\t{ref}\t{alt}\t99\tPASS\t{info}\tGT\t0/1\n"
                   for p, i, ref, alt, info in records)
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.vcf"
        src.write_text(HEADER + "\n" + body)
        r = subprocess.run(["awk", "-f", str(AWK), str(src)],
                           capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"awk failed: {r.stderr.strip()[:300]}")
    rows = [l.split("\t") for l in r.stdout.splitlines() if not l.startswith("#")]
    return rows, r.stderr


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("a sequence-resolved ALT becomes symbolic")
    # The shape from the real failure: Manta writes out the sequence it assembled.
    rows, _ = run([(1000, "r1", "TTTTTTCTTTCTT", "T", "SVTYPE=DEL;END=1012")])
    check("ALT is replaced by <SVTYPE>", rows[0][4] == "<DEL>", f"ALT={rows[0][4]}")
    check("REF is cut to the single anchoring base", rows[0][3] == "T",
          "a long REF beside a symbolic ALT states two different lengths for one record")
    check("INFO is untouched", rows[0][7] == "SVTYPE=DEL;END=1012", f"INFO={rows[0][7]}")

    rows, _ = run([(1000, "r1", "T", "TCTTTCTTTCTTTCTT", "SVTYPE=INS;END=1000")])
    check("an insertion's assembled sequence is replaced too", rows[0][4] == "<INS>",
          f"ALT={rows[0][4]}")

    print("everything SVAnnotate already accepts is left alone")
    rows, _ = run([(1000, "r1", "N", "<DEL>", "SVTYPE=DEL;END=2000"),
                   (2000, "r2", "N", "<DUP:TANDEM>", "SVTYPE=DUP;END=3000"),
                   (3000, "r3", "N", "N[chr2:5000[", "SVTYPE=BND"),
                   (4000, "r4", "N", ".", "SVTYPE=DEL;END=4100")])
    check("a symbolic ALT is unchanged", rows[0][4] == "<DEL>", f"ALT={rows[0][4]}")
    check("a subtyped symbolic ALT keeps its subtype",
          rows[1][4] == "<DUP:TANDEM>",
          f"rewriting it from SVTYPE would flatten it to <DUP>: {rows[1][4]}")
    check("a breakend is unchanged", rows[2][4] == "N[chr2:5000[", f"ALT={rows[2][4]}")
    check("a missing ALT is unchanged", rows[3][4] == ".", f"ALT={rows[3][4]}")
    check("REF is not cut on records that were left alone",
          [r[3] for r in rows] == ["N", "N", "N", "N"])

    print("a multi-allelic record is never collapsed into one symbol")
    # The silent-wrong case this guard exists for. A record carries ONE SVTYPE, so writing
    # <DEL> over the whole column would turn two distinct alleles into one -- and would do it
    # without erroring, which is the failure class this pipeline keeps trying to avoid.
    rows, err = run([(1000, "r1", "A", "TTTTTTCTT,<DEL>", "SVTYPE=DEL;END=1012")])
    check("the ALT column is left exactly as it was",
          rows[0][4] == "TTTTTTCTT,<DEL>", f"ALT={rows[0][4]}")
    check("REF is left alone too", rows[0][3] == "A", f"REF={rows[0][3]}")
    check("and stderr says SVAnnotate will reject it",
          "multi-allelic" in err and "WARNING" in err, f"stderr={err.strip()[:200]}")

    rows, err = run([(1000, "r1", "A", "A,TTTTTTCTT", "SVTYPE=DEL;END=1012")])
    check("a literal allele in second position is caught too",
          rows[0][4] == "A,TTTTTTCTT" and "multi-allelic" in err,
          "scanning only the first allele would miss this one")

    rows, err = run([(1000, "r1", "N", "<DEL>,<DUP>", "SVTYPE=DEL;END=1012")])
    check("an all-symbolic multi-allelic record is unchanged and unwarned",
          rows[0][4] == "<DEL>,<DUP>" and err.strip() == "",
          f"SVAnnotate accepts these already: stderr={err.strip()[:200]}")

    print("a literal ALT with no SVTYPE is reported, not guessed at")
    rows, err = run([(1000, "r1", "TTTTTT", "T", "END=1005")])
    check("the record passes through untouched", rows[0][4] == "T" and rows[0][3] == "TTTTTT",
          "inferring DEL from REF/ALT lengths would be inventing a call")
    check("and stderr says SVAnnotate will still reject it",
          "no SVTYPE" in err and "WARNING" in err, f"stderr={err.strip()[:200]}")

    print("the rewrite is counted, so a silent no-op is visible")
    _, err = run([(1000, "r1", "TTTTTT", "T", "SVTYPE=DEL;END=1005")])
    check("stderr reports how many were rewritten", "rewrote 1" in err,
          f"stderr={err.strip()[:200]}")
    _, err = run([(1000, "r1", "N", "<DEL>", "SVTYPE=DEL;END=2000")])
    check("and says nothing when there was nothing to do", err.strip() == "",
          f"stderr={err.strip()[:200]}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("symbolic ALT rewrite holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
