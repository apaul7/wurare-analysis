#!/usr/bin/env python3
"""annotate_snps' --alignments path: somalier extracts from CRAM/BAMs, not the VCF.

    python3 tests/test_somalier_alignments.py

Runs the real DAG with test_skip_params' fake tools (imported, not copied -- the pattern
test_summary_report.py established). What is under test is the WIRING: with a sheet, the
per-sample somalier_extract runs and somalier_extract_vcf must not; everything downstream
(relate, check_somalier, the summary report's qc rows) is fed the same .somalier files and
must be unaware of the switch. The fake somalier's `extract` branch already writes
SAMP_M.somalier and SAMP_F.somalier into -d, which satisfies somalier_extract's per-sample
"[ -f <sample>.somalier ]" guard as long as the sheet uses those names.

Real-tool depth semantics for alignment extraction are pinned by tests/test_ploidy_e2e.py
(annotate_svs, real somalier in Docker) and are not duplicated here. The VCF path's
continued operation is pinned by test_skip_params' own QC assertions.

Also covered: the startup guards -- a sheet without the alignment reference, without its
index, or without --somalier_sites must be rejected at startup by name, and a sheet
missing a column must name the column.
"""

import shutil
import sys
import tempfile
from pathlib import Path

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def main():
    if not shutil.which("nextflow"):
        print("SKIP: nextflow not on PATH")
        return 0
    import test_skip_params as tsp

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        binf = tsp.build_env(tmp)

        # Alignment fixtures. Content-free: the fake somalier never opens them, and the
        # pipeline only checkIfExists's them on the head node.
        for name in ("SAMP_M.bam", "SAMP_M.bam.bai", "SAMP_F.bam", "SAMP_F.bam.bai"):
            (tmp / name).touch()
        (tmp / "alignments.csv").write_text(
            "sample,alignment,alignment_index\n"
            f"SAMP_M,{tmp}/SAMP_M.bam,{tmp}/SAMP_M.bam.bai\n"
            f"SAMP_F,{tmp}/SAMP_F.bam,{tmp}/SAMP_F.bam.bai\n")

        aln_args = ["--alignments", "alignments.csv",
                    "--alignment_reference", "ref.fa",
                    "--alignment_reference_index", "ref.fa.fai"]

        # --- the sheet switches the extract source -----------------------------------
        # --skip_spliceai_squirls so the run completes locally (see test_skip_params).
        print("with a sheet: somalier_extract runs per sample, the VCF extract does not")
        trace, _named, out = tsp.run_pipeline(
            tmp, binf,
            ["--skip_spliceai_squirls",
             "--somalier_sites", "sites.vcf.gz", "--ped", "cohort.ped"] + aln_args, "aln")
        check("the run completes",
              "[FAILED]" not in out and "ERROR" not in out, f"...{out.strip()[-500:]}")
        check("somalier_extract ran", "somalier_extract" in trace, str(sorted(trace)))
        check("somalier_extract_vcf did NOT run",
              "somalier_extract_vcf" not in trace, str(sorted(trace)))
        for proc in ("somalier_relate", "check_somalier"):
            check(f"{proc} still ran", proc in trace, str(sorted(trace)))

        qc_dir = tmp / "results_aln" / "04_qc"
        ploidy = qc_dir / "ploidy.tsv"
        check("04_qc/ploidy.tsv is published", ploidy.is_file(),
              str(sorted(p.name for p in qc_dir.glob("*")) if qc_dir.is_dir()
                  else "no 04_qc"))
        if ploidy.is_file():
            rows = {r[0]: r[1:] for r in
                    (l.split("\t") for l in ploidy.read_text().splitlines()
                     if not l.startswith("#") and l.strip())}
            check("the hemizygous sample is still called 1 X / 1 Y",
                  rows.get("SAMP_M", [])[:2] == ["1", "1"], str(rows))
            check("the diploid sample is still called 2 X / 0 Y",
                  rows.get("SAMP_F", [])[:2] == ["2", "0"], str(rows))

        summary = tmp / "results_aln" / "05_report" / "TEST_wgs_20260627.summary.tsv"
        check("05_report summary is published", summary.is_file(), str(summary))
        if summary.is_file():
            srows = {}
            for line in summary.read_text().splitlines()[1:]:
                section, metric, sample, value = line.split("\t")
                srows[(section, metric, sample)] = value
            check("the report's inferred sex comes through the alignment path",
                  srows.get(("qc", "inferred_sex", "SAMP_M")) == "male"
                  and srows.get(("qc", "inferred_sex", "SAMP_F")) == "female", str(srows))

        # --- startup guards ----------------------------------------------------------
        print("a sheet without its reference, index or sites file is rejected at startup")
        for missing, args, tag in [
            ("alignment_reference",
             ["--somalier_sites", "sites.vcf.gz",
              "--alignments", "alignments.csv",
              "--alignment_reference_index", "ref.fa.fai"], "noref"),
            ("alignment_reference_index",
             ["--somalier_sites", "sites.vcf.gz",
              "--alignments", "alignments.csv",
              "--alignment_reference", "ref.fa"], "noidx"),
            ("somalier_sites", aln_args, "nosites"),
        ]:
            trace, _named, out = tsp.run_pipeline(tmp, binf, args, tag)
            check(f"missing {missing} is named in the error",
                  missing in out, f"...{out.strip()[-300:]}")
            check(f"and nothing was submitted ({tag})", not trace, str(sorted(trace)))

        print("a sheet missing a column is rejected naming the column")
        (tmp / "bad.csv").write_text(
            f"sample,alignment\nSAMP_M,{tmp}/SAMP_M.bam\n")
        trace, _named, out = tsp.run_pipeline(
            tmp, binf,
            ["--somalier_sites", "sites.vcf.gz",
             "--alignments", "bad.csv",
             "--alignment_reference", "ref.fa",
             "--alignment_reference_index", "ref.fa.fai"], "badcol")
        check("the error names alignment_index", "alignment_index" in out,
              f"...{out.strip()[-300:]}")
        check("and nothing was submitted (badcol)", not trace, str(sorted(trace)))

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("the alignments path holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
