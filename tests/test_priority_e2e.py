#!/usr/bin/env python3
"""End-to-end regression test for --priority order at both merge axes.

    python tests/test_priority_e2e.py

Like test_grouping_e2e.py this runs the real pipeline, because the thing under test is the
order in which merge_svs hands files to SVDB and that lives in Nextflow channel wiring. It
skips cleanly when nextflow, bcftools or docker are unavailable.

Both assertions are grounded in what the design spike MEASURED: where two inputs
genotype the same sample differently, the --priority winner's
genotype is the one on the merged record. That makes the genotype a direct readout of the
priority order, which is otherwise invisible in the output.

  axis A   manta first, then the rest alphabetically (decided 2026-07-31). Two
           callers over the SAME sample set meet here, and manta's call is the one to keep.

  axis B   joint inputs first. A sample present in both a joint VCF and its own
           single-sample calls must keep the JOINT genotype -- that is the entire reason
           the design refuses to split multi-sample VCFs. This is the case that was wrong: axis B
           sorted on the label alone, so priority went to whichever sample_set sorted first
           and the joint VCF won or lost by luck of the alphabet.

The genotypes are chosen so that a wrong order produces a DIFFERENT value rather than an
absent one -- 0/1 vs 1/1, both valid, neither missing. A test that only checks "a genotype
is present" passes under either order.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
MAIN = ROOT / "nf" / "annotate_svs" / "main.nf"

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


def write_vcf(tmp, name, genotypes, end=2000, svlen=-1000):
    """genotypes: list of (sample, GT). One DEL record at chr1:1000."""
    plain = tmp / f"{name}.vcf"
    samples = "\t".join(s for s, _ in genotypes)
    calls = "\t".join(gt for _, gt in genotypes)
    plain.write_text(
        HEADER
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n"
        + f"chr1\t1000\t.\tN\t<DEL>\t99\tPASS\t"
          f"SVTYPE=DEL;END={end};SVLEN={svlen}\tGT\t{calls}\n")
    gz = tmp / f"{name}.vcf.gz"
    gz.write_bytes(subprocess.run(["bgzip", "-c", str(plain)], capture_output=True).stdout)
    run(["tabix", "-p", "vcf", str(gz)])
    return gz


def run_pipeline(tmp, sheet_rows, samples, env):
    """sheet_rows: list of (sample_set, caller, joint, vcf). Returns {sample: GT} or None."""
    sheet = tmp / "vcfs.csv"
    sheet.write_text("sample_set,caller,joint,vcf,tbi\n" + "".join(
        f"{label},{caller},{str(joint).lower()},{gz},{gz}.tbi\n"
        for label, caller, joint, gz in sheet_rows))
    ped = tmp / "cohort.ped"
    ped.write_text("".join(f"FAM_{s}\t{s}\t0\t0\t1\t0\n" for s in samples))

    r = run(["nextflow", "-q", "run", str(MAIN), "-with-docker",
             "--vcfs", str(sheet), "--ped", str(ped),
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

    q = run(["bcftools", "query", "-f", "[%SAMPLE=%GT\\n]", str(cohort)])
    lines = [l for l in q.stdout.splitlines() if l]
    # One line per sample per record, so more lines than samples means more than one record
    # survived -- and on unmerged records the priority order is not observable at all.
    check("the two inputs collapsed to ONE record",
          len(lines) == len(set(l.split("=")[0] for l in lines)),
          f"got {lines}")
    return dict(l.split("=", 1) for l in lines)


def main():
    for tool in ("nextflow", "bcftools", "bgzip", "tabix", "docker"):
        if not shutil.which(tool):
            print(f"SKIP: {tool} not on PATH "
                  "(needs the annotate-svs conda env and a running Docker)")
            return 0
    if run(["docker", "info"]).returncode != 0:
        print("SKIP: docker daemon not running")
        return 0

    env = dict(os.environ, NXF_SYNTAX_PARSER="v2")

    print("axis A -- manta outranks an alphabetically earlier caller")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # Same sample set, so these meet at axis A. 700/1000 = 0.70 reciprocal overlap,
        # above --overlap_axis_a (0.6). delly sorts first alphabetically, which is exactly
        # what the manta-first rule has to override.
        manta = write_vcf(tmp, "manta_s1", [("SAMP1", "0/1")], end=2000, svlen=-1000)
        delly = write_vcf(tmp, "delly_s1", [("SAMP1", "1/1")], end=1700, svlen=-700)
        gts = run_pipeline(tmp, [("SET", "manta", False, manta),
                                 ("SET", "delly", False, delly)], ["SAMP1"], env)
        if gts is not None:
            check("the merged record keeps manta's genotype",
                  gts.get("SAMP1") == "0/1",
                  f"SAMP1={gts.get('SAMP1')} -- 1/1 means delly took priority, which is "
                  "what plain alphabetical ordering does")

    print("axis B -- a joint VCF outranks the same sample's single-sample calls")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # Different sample sets, so neither meets the other at axis A and both arrive at
        # axis B. Identical coordinates, so they certainly cluster at --overlap_axis_b.
        # The joint set is labelled ZZZ and the single AAA: under the old label-only sort
        # the single-sample VCF took priority, which is the defect this pins.
        joint = write_vcf(tmp, "joint_set", [("SAMP1", "1/1"), ("SAMP2", "0/1")])
        single = write_vcf(tmp, "manta_s1", [("SAMP1", "0/1")])
        gts = run_pipeline(tmp, [("ZZZ_JOINT", "gatk", True, joint),
                                 ("AAA_SINGLE", "manta", False, single)],
                           ["SAMP1", "SAMP2"], env)
        if gts is not None:
            check("the shared sample keeps the JOINT genotype",
                  gts.get("SAMP1") == "1/1",
                  f"SAMP1={gts.get('SAMP1')} -- 0/1 is the single-sample VCF's call, which "
                  "means priority went by label and the design's reason for keeping joint VCFs "
                  "unsplit was silently discarded")
            check("the joint-only sample survives the merge",
                  gts.get("SAMP2") == "0/1", f"SAMP2={gts.get('SAMP2')}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("priority order holds at both axes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
