#!/usr/bin/env python3
"""End-to-end regression test for the axis-A grouping key.

Unlike the other suites this one runs the actual pipeline, because the defect it guards
against lives in Nextflow channel wiring and cannot be reproduced any other way. It skips
cleanly when nextflow, bcftools or docker are unavailable.

    python tests/test_grouping_e2e.py

The bug it caught: merge_svs grouped axis A on the samplesheet's `sample_set` column -- free
text a human types -- while the design, the module's
own comment, and the README all say the key is the SET of sample IDs read from the VCF
header. prepare_svs computed the correct key, printed it at preflight, and threw it away.

The consequence is subtler than it first looks, and worth stating precisely because the
obvious test does NOT catch it. When axis A never groups two callers, their records still
reach axis B and still cluster there, and promote_caller_support.awk in cohort mode still
recovers CALLER_SUPP from the per-input <tag>_INFO blobs. So on a well-overlapping call the
output is identical either way.

What actually differs is the THRESHOLD applied. Axis A is cross-caller over the same samples
and runs at --overlap_axis_a (0.6); axis B is cross-sample assembly at --overlap_axis_b
(0.8). A pair of calls overlapping between those two values merges under the correct
grouping and does not under the buggy one. This fixture is built at ~0.70 reciprocal overlap
precisely to sit in that window -- at 0.99, as the committed spike fixtures do, both
groupings agree and the bug is invisible.

Two scenarios, and they fail in opposite directions -- which is why both are needed. A test
that only checks "did these merge" is satisfied by a grouping that merges everything.

  same sample, different labels   must merge   (grouping on the label wrongly splits them)
  different samples, same label   must NOT     (grouping on the label wrongly joins them)

The second is the worse bug: merging across samples at cross-caller thresholds makes a
one-caller/two-sample shared event indistinguishable from a two-caller/one-sample duplicate,
which is exactly the confusion the two-axis design exists to prevent. It is a wrong answer
rather than a missing one.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
MAIN = ROOT / "nf" / "annotate_svs" / "main.nf"

# Deliberately 700/1000 = 0.70 reciprocal overlap: between --overlap_axis_a (0.6) and
# --overlap_axis_b (0.8). That window is the ONLY place the grouping key is observable --
# outside it both groupings agree and the bug is invisible, which is how the first version
# of this test passed against the bug.
HEADER = ("##fileformat=VCFv4.2\n"
          "##contig=<ID=chr1,length=248956422>\n"
          "##ALT=<ID=DEL,Description=\"Deletion\">\n"
          "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n"
          "##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n"
          "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"l\">\n"
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n")

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def write_vcf(tmp, name, sample, end, svlen):
    plain = tmp / f"{name}.vcf"
    plain.write_text(
        HEADER
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        + f"chr1\t1000\t.\tN\t<DEL>\t99\tPASS\t"
          f"SVTYPE=DEL;END={end};SVLEN={svlen}\tGT\t0/1\n")
    gz = tmp / f"{name}.vcf.gz"
    gz.write_bytes(subprocess.run(["bgzip", "-c", str(plain)], capture_output=True).stdout)
    run(["tabix", "-p", "vcf", str(gz)])
    return gz


def write_ped(tmp, samples):
    """The PED is a required input, so every e2e run needs one. Sex 1 (male) throughout --
    these fixtures live on chr1, where sex does not enter the depth logic."""
    ped = tmp / "cohort.ped"
    ped.write_text("".join(f"FAM_{s}\t{s}\t0\t0\t1\t0\n" for s in samples))
    return ped


def run_pipeline(tmp, sheet_rows, env, samples):
    """sheet_rows: list of (sample_set, caller, vcf_path). Returns cohort rows or None."""
    sheet = tmp / "vcfs.csv"
    sheet.write_text("sample_set,caller,joint,vcf,tbi\n" + "".join(
        f"{label},{caller},false,{gz},{gz}.tbi\n" for label, caller, gz in sheet_rows))
    r = run(["nextflow", "-q", "run", str(MAIN), "-with-docker",
             "--vcfs", str(sheet), "--ped", str(write_ped(tmp, samples)),
             "--outdir", str(tmp / "results")],
            cwd=str(ROOT), env=env)
    errors = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith("[ERROR]")]
    check("the pipeline runs clean", not errors, "; ".join(errors[:2]))
    if errors:
        return None
    cohort = tmp / "results" / "05_filter" / "cohort.tagged.vcf.gz"
    if not cohort.is_file():
        check("a cohort VCF was produced", False, str(cohort))
        return None
    q = run(["bcftools", "query", "-f",
             "%CHROM:%POS\t%INFO/CALLER_SUPP\t%INFO/NCALLER\t[%GT ]\n", str(cohort)])
    # A LIST, not a dict keyed on position: these records share a position, and a
    # position-keyed dict silently collapses them so the count assertion passes vacuously.
    # (It did, on the first version of this test.)
    return [tuple(line.split("\t")) for line in q.stdout.splitlines()]


def main():
    for tool in ("nextflow", "bcftools", "bgzip", "tabix", "docker"):
        if not shutil.which(tool):
            print(f"SKIP: {tool} not on PATH "
                  "(needs the annotate-svs conda env and a running Docker)")
            return 0
    if run(["docker", "info"]).returncode != 0:
        print("SKIP: docker daemon not running")
        return 0

    env = dict(os.environ, NXF_SYNTAX_PARSER="v2", NXF_ANSI_LOG="false")

    print("same sample, different sheet labels -- must MERGE at axis A")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        a = write_vcf(tmp, "manta_one", "SAMP1", 2000, -1000)
        b = write_vcf(tmp, "delly_one", "SAMP1", 1700, -700)
        rows = run_pipeline(tmp, [("LABEL_A", "manta", a), ("LABEL_B", "delly", b)], env,
                            ["SAMP1"])
        if rows is not None:
            check("the two callers collapse to ONE record at 0.70 overlap",
                  len(rows) == 1,
                  f"got {len(rows)} ({rows}) -- two means axis A never grouped them and "
                  "they fell through to axis B's tighter threshold")
            if len(rows) == 1:
                _pos, supp, n, _gt = rows[0]
                check("and it carries both callers",
                      sorted(supp.split(",")) == ["delly", "manta"], f"CALLER_SUPP={supp}")
                check("NCALLER counts both", n == "2", f"NCALLER={n}")

    print("different samples, same sheet label -- must NOT merge across samples")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        a = write_vcf(tmp, "manta_s1", "SAMP1", 2000, -1000)
        b = write_vcf(tmp, "delly_s2", "SAMP2", 1700, -700)
        # Same label on two DIFFERENT samples. Grouping on the label puts them in one axis-A
        # group and merges them at 0.6; grouping on the sample set keeps them apart, and
        # axis B's 0.8 then correctly refuses to merge a 0.70 overlap.
        rows = run_pipeline(tmp, [("SHARED", "manta", a), ("SHARED", "delly", b)], env,
                            ["SAMP1", "SAMP2"])
        if rows is not None:
            check("two samples' distinct calls stay TWO records",
                  len(rows) == 2,
                  f"got {len(rows)} ({rows}) -- one means axis A merged across samples at "
                  "cross-caller thresholds, which makes a one-caller/two-sample shared "
                  "event indistinguishable from a two-caller/one-sample duplicate")
            check("neither record claims corroboration it does not have",
                  all(n == "1" for _p, _s, n, _g in rows),
                  f"rows={rows} -- NCALLER=2 here would be a fabricated agreement between "
                  "two different samples")
            check("both samples appear in the cohort",
                  len(rows) == 2 and {r[1] for r in rows} == {"manta", "delly"},
                  f"rows={rows}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("axis-A grouping holds in both directions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
