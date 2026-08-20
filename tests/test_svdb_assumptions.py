#!/usr/bin/env python3
"""Regression tests for the SVDB behaviours the pipeline design depends on.

These are characterization tests against the *tool*, not against pipeline code
(there is no pipeline yet). Each one pins an assumption the design spike established by
experiment. If an SVDB upgrade changes any of them, the concern named in
the failure message needs revisiting before the upgrade ships.

Some assert behaviour that is *wrong* but currently true -- the duplicate-INFO-
key defect, and ALGORITHMS under-reporting. Those are marked KNOWN-DEFECT: when
they fail, upstream has fixed it and the pipeline stage that works around it
can go away.

    python tests/test_svdb_assumptions.py

Needs `svdb` and `bcftools` on PATH (conda env `annotate-svs`, environment.yml).
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures" / "svdb_spike"
INFO_COL, FORMAT_COL = 7, 8

failures = []


def check(name, section, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  [{section}]  {detail}")
        failures.append(name)


def run(cmd):
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True)


def merge(out, inputs):
    """inputs: list of (path, tag), highest priority first."""
    cmd = ["svdb", "--merge", "--vcf"] + [f"{p}:{t}" for p, t in inputs]
    cmd += ["--priority", ",".join(t for _, t in inputs)]
    res = run(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"svdb merge failed: {res.stderr.strip()[:300]}")
    out.write_text(res.stdout)
    return out


def records(vcf):
    return [ln.split("\t") for ln in vcf.read_text().splitlines() if not ln.startswith("#")]


def sample_names(vcf):
    line = next(ln for ln in vcf.read_text().splitlines() if ln.startswith("#CHROM"))
    return line.split("\t")[FORMAT_COL + 1:]


def info_keys(record):
    return [kv.split("=")[0] for kv in record[INFO_COL].split(";") if "=" in kv]


def main():
    for tool in ("svdb", "bcftools"):
        if not shutil.which(tool):
            # Skip, not fail, matching every other test here: a missing conda env says
            # nothing about the code. This exited non-zero until tests/run.sh made the
            # difference visible as the suite's one red line.
            print(f"SKIP: {tool} not on PATH -- activate the annotate-svs conda env")
            return 0
    missing = [v for v in ("manta_SAMP1.vcf", "delly_SAMP1.vcf", "joint.vcf", "samp3.vcf")
               if not (FIXTURES / v).is_file()]
    if missing:
        sys.exit(f"fixtures missing from {FIXTURES}: {', '.join(missing)}")

    with tempfile.TemporaryDirectory() as td:
        tmp, f = Path(td), FIXTURES

        # --- axis A: two callers, one sample -------------------------------
        print("axis A (cross-caller, one sample)")
        axis_a = merge(tmp / "axisA.vcf",
                       [(f / "manta_SAMP1.vcf", "manta"), (f / "delly_SAMP1.vcf", "delly")])
        recs_a = {(r[0], r[1]): r for r in records(axis_a)}

        check("shared DEL collapses to one record", "two-axis merge",
              len(recs_a) == 4,
              f"expected 4 (1 merged + 3 singletons), got {len(recs_a)}")

        shared_a = recs_a[("chr1", "1000")][INFO_COL]
        check("merged record names both callers in svdb_origin", "provenance promotion",
              "svdb_origin=manta|delly" in shared_a,
              "svdb_origin is the CALLER_SUPP source; if renamed, the promotion stage changes")
        check("merged record counts both callers in FOUNDBY", "provenance promotion",
              "FOUNDBY=2" in shared_a, "FOUNDBY is the NCALLER source")
        check("KNOWN-DEFECT: ALGORITHMS keeps only the priority caller", "ALGORITHMS re-derivation",
              "ALGORITHMS=manta;" in shared_a + ";",
              "if this fails, upstream fixed it and the ALGORITHMS re-derivation is redundant")

        # SUPP_VEC bit order is alphabetical by tag, NOT --vcf order. Inputs were
        # given manta,delly; a delly-only record must be '10' because delly<manta.
        check("SUPP_VEC bits are alphabetical by tag, not input order", "spike finding; soft filters",
              "SUPP_VEC=10" in recs_a[("chr3", "1000")][INFO_COL],
              "never read SUPP_VEC positionally against the --vcf list")

        # --- axis B: cohort assembly ---------------------------------------
        print("axis B (cohort assembly)")
        axis_b = merge(tmp / "axisB.vcf",
                       [(f / "joint.vcf", "jointcaller"), (axis_a, "samp1set"),
                        (f / "samp3.vcf", "samp3set")])
        cols = sample_names(axis_b)
        recs_b = {(r[0], r[1]): r for r in records(axis_b)}
        gt = {k: dict(zip(cols, r[FORMAT_COL + 1:])) for k, r in recs_b.items()}

        check("sample columns are unioned across inputs", "cohort matrix",
              cols == ["SAMP1", "SAMP2", "SAMP3"], f"got {cols}")

        # THE claim the whole plan rests on: SAMP3's input loses the cluster to
        # the joint VCF, and its genotype must survive anyway.
        shared_b = gt[("chr1", "1000")]
        check("genotype survives for a sample whose input LOST the cluster", "spike finding",
              shared_b["SAMP3"] == "0/1",
              f"SAMP3={shared_b['SAMP3']} -- if './.', axis B is not a cohort matrix "
              "and the two-axis recommendation is void")
        check("priority winner supplies the genotype for a shared sample", "genotype priority",
              shared_b["SAMP1"] == "0/1" and shared_b["SAMP2"] == "1/1", f"got {shared_b}")

        joint_only = gt[("chr4", "3000")]
        check("0/0 stays distinct from ./.", "genotype priority",
              joint_only["SAMP2"] == "0/0" and joint_only["SAMP3"] == "./.",
              f"got {joint_only} -- this distinction is what the joint-VCF handling protects")

        check("<tag>_SAMPLE records which input supplied each genotype", "spike finding; provenance",
              "jointcaller_SAMPLE=joint_1|SAMP1|GT:0/1|SAMP2|GT:1/1"
              in recs_b[("chr1", "1000")][INFO_COL],
              "if this format changes, the provenance mitigation needs rework")

        # --- the defect ----------------------------------------------------
        print("stacked-merge INFO key collision")
        keys = info_keys(recs_b[("chr1", "1000")])
        dup = sorted(k for k in ("svdb_origin", "FOUNDBY", "SUPP_VEC", "set", "VARID")
                     if keys.count(k) > 1)
        check("KNOWN-DEFECT: axis B duplicates axis A's INFO keys", "spike finding; provenance promotion",
              dup == ["FOUNDBY", "SUPP_VEC", "VARID", "set", "svdb_origin"],
              f"duplicated: {dup or 'none'} -- if none, upstream fixed it and "
              "the strip stage can be dropped")

        q = run(["bcftools", "query", "-f", "%INFO/FOUNDBY %INFO/SUPP_VEC\n", axis_b])
        first = q.stdout.splitlines()[0].strip() if q.stdout else ""
        check("KNOWN-DEFECT: bcftools returns the stale axis-A value", "spike finding",
              first == "2 11",
              f"got {first!r}, expected '2 11' (axis-B truth is '3 111')")

        # --- invocation contract -------------------------------------------
        print("invocation contract")
        untagged = run(["svdb", "--merge", "--vcf", f / "joint.vcf", f / "samp3.vcf",
                        "--priority", "jointcaller,samp3set"])
        check("--priority on untagged inputs is a hard error, not a silent no-op", "two-axis merge",
              untagged.returncode != 0 or "mismatch" in untagged.stdout + untagged.stderr,
              "if this becomes silent, the merge stage needs a defensive assertion again")

        # Priority order decides the genotype, not just the record shape.
        conflict = tmp / "conflict.vcf"
        conflict.write_text(axis_a.read_text().replace("GT\t0/1", "GT\t1/1"))
        joint_first = merge(tmp / "p1.vcf",
                            [(f / "joint.vcf", "jointcaller"), (conflict, "samp1set")])
        single_first = merge(tmp / "p2.vcf",
                             [(conflict, "samp1set"), (f / "joint.vcf", "jointcaller")])
        g1 = next(r for r in records(joint_first) if r[1] == "1000")[FORMAT_COL + 1]
        g2 = next(r for r in records(single_first) if r[1] == "1000")[FORMAT_COL + 1]
        check("--priority order changes which genotype survives", "provenance",
              g1 == "0/1" and g2 == "1/1",
              f"joint-first={g1}, single-first={g2} -- if equal, priority stopped "
              "mattering and the joint-VCF protection needs another mechanism")

        # --- thresholds -----------------------------------------------------
        print("threshold defaults")
        helptext = run(["svdb", "--merge", "--help"]).stdout
        check("--overlap default is still 0.95", "two-axis merge; provenance", "0.95" in helptext)
        check("--bnd_distance default is still 2000", "two-axis merge", "2000" in helptext)
        # --ins_distance moved 50 -> 25 between 2.8.1 and 2.12.0. No assertion on
        # the value; the point is that it drifts, so the pipeline must set it.
        print("  note  --ins_distance drifts across releases (50 in 2.8.1, 25 in "
              "2.12.0)\n        -- pin the version, set all three thresholds explicitly")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        print("A failure means a plan assumption no longer holds. Read the cited "
              "section before changing anything.")
        return 1
    print("all assumptions hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
