#!/usr/bin/env python3
"""Regression tests for the AnnotSV coverage report.

Covers nf/annotate_svs/assets/check_annotsv_coverage.awk directly -- no Nextflow, no
container, and notably no AnnotSV, which needs a multi-gigabyte annotation bundle.

The design originally asked for "AnnotSV TSV row count matches the cohort VCF
record count". That check fails on a correct run: AnnotSV's default -annotationMode is
`both` (full AND split rows), and even in `full` mode it legitimately drops records below
its -SVminSize and records it cannot type. What matters is WHICH records went missing, so
the report lists them by ID instead of asserting a number.

    python tests/test_annotsv_coverage.py

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

AWK = (Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
       / "check_annotsv_coverage.awk")

VCF_HEADER = ("##fileformat=VCFv4.2\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1")

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def coverage(vcf_ids, tsv_ids, tsv_columns=("AnnotSV_ID", "SV_type", "ID", "ACMG_class")):
    vcf = VCF_HEADER + "\n" + "\n".join(
        f"chr1\t{1000 + i}\t{rid}\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL\tGT\t0/1"
        for i, rid in enumerate(vcf_ids))
    id_at = tsv_columns.index("ID") if "ID" in tsv_columns else None
    rows = []
    for rid in tsv_ids:
        cells = [f"x{c}" for c in range(len(tsv_columns))]
        if id_at is not None:
            cells[id_at] = rid
        rows.append("\t".join(cells))
    tsv = "\t".join(tsv_columns) + "\n" + "\n".join(rows)
    with tempfile.TemporaryDirectory() as td:
        v = Path(td) / "cohort.vcf"; v.write_text(vcf + "\n")
        t = Path(td) / "annotsv.tsv"; t.write_text(tsv + "\n")
        r = subprocess.run(["awk", "-f", str(AWK), str(v), str(t)],
                           capture_output=True, text=True)
    return r


def value(out, label):
    for line in out.splitlines():
        if line.startswith(label + "\t"):
            return line.split("\t")[1]
    return None


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("full coverage")
    r = coverage(["a", "b", "c"], ["a", "b", "c"])
    check("exits 0", r.returncode == 0, r.stderr.strip()[:120])
    check("counts the cohort records", value(r.stdout, "cohort VCF records") == "3")
    check("counts the TSV rows", value(r.stdout, "AnnotSV TSV rows") == "3")
    check("reports nothing absent",
          value(r.stdout, "cohort records absent from TSV") == "0", r.stdout)

    print("records AnnotSV dropped are named, not just counted")
    r = coverage(["a", "b", "c"], ["a", "c"])
    check("counts the absentees", value(r.stdout, "cohort records absent from TSV") == "1")
    check("names the absent ID", "\nb" in r.stdout,
          "a count alone does not tell the reader which record vanished")
    check("still exits 0", r.returncode == 0,
          "a dropped record is information, not a pipeline failure")

    print("split-mode style duplicate rows do not inflate the match")
    # -annotationMode both emits one row per SV x gene; the same ID appears many times.
    r = coverage(["a", "b"], ["a", "a", "a", "b"])
    check("TSV row count reflects the duplicates",
          value(r.stdout, "AnnotSV TSV rows") == "4")
    check("distinct cohort records matched is still 2",
          value(r.stdout, "cohort records present in TSV") == "2",
          "counting rows instead of distinct IDs is exactly the check that was wrong")
    check("nothing reported absent",
          value(r.stdout, "cohort records absent from TSV") == "0")

    print("a mismatched pairing is surfaced")
    r = coverage(["a", "b"], ["a", "zzz"])
    check("an unexpected TSV ID is reported", "zzz" in r.stdout,
          "means the TSV and VCF are not from the same run")
    check("and the missing cohort record is too", "\nb" in r.stdout)

    print("the ID column is found by name, not position")
    r = coverage(["a"], ["a"], tsv_columns=("ID", "SV_type"))
    check("ID first still works", value(r.stdout, "cohort records present in TSV") == "1")
    r = coverage(["a"], ["a"], tsv_columns=("AnnotSV_ID", "SV_type", "Gene_name", "ID"))
    check("ID last still works", value(r.stdout, "cohort records present in TSV") == "1",
          "AnnotSV column order differs between releases")

    print("a TSV without an ID column is an error, not a silent pass")
    r = coverage(["a"], ["a"], tsv_columns=("AnnotSV_ID", "SV_type"))
    check("exits non-zero", r.returncode != 0,
          "without an ID column the comparison never happened, so reporting success "
          "would be worse than failing")
    check("says why", "no 'ID' column" in r.stderr, r.stderr.strip()[:120])

    print("empty inputs")
    r = coverage([], [])
    check("no records anywhere still exits 0", r.returncode == 0)
    check("reports zero cohort records", value(r.stdout, "cohort VCF records") == "0",
          r.stdout)

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("coverage reporting holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
