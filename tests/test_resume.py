#!/usr/bin/env python3
"""End-to-end regression test for `-resume`.

    python tests/test_resume.py

The defect it guards against: -resume cached NOTHING, on every run, for months. A pipeline
this slow is expected to be restartable, and it silently was not -- a run that died in the
cohort merge re-did Phase 1 from the top.

The cause was records in process `input:` blocks. Nextflow has no value-based hash for a
`record`, so its task hasher falls back to Object.toString(), which embeds the JVM identity
hash:

    9cb09c67... [Main.Scalars] Main.Scalars@217bf527

That value is a fresh allocation address every run, so the task hash changed on every run and
could never match a cache entry. Records are fine in `output:` blocks and fine in channel
plumbing -- only `input:` is hashed. Processes therefore take label/vcf/tbi/joint and the
workflows destructure at the call site.

WHY THE FIXTURE IS BUILT ONCE AND BOTH RUNS SHARE IT, which is the whole trick of this test.
Nextflow's default cache mode hashes each input file by path, size and LAST-MODIFIED. The
other e2e suites build their fixture in a fresh `tempfile.TemporaryDirectory()` per
invocation, so running one of them twice produces new paths and new mtimes and therefore
`cached=0` -- with or without the bug. Written that way this test would fail forever, for a
reason that has nothing to do with what it is testing. So the fixture is written once, and
both runs read those same files.

Both runs also launch from inside the fixture directory, so `.nextflow/` (which holds the
cache DB `-resume` reads) and `work/` land there and are cleaned up with it. Launching from
the repo root would leave the run's cache sitting in the working tree.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from test_skip_params import trace_stats

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


def write_vcf(tmp, name, samples):
    plain = tmp / f"{name}.vcf"
    plain.write_text(
        HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
        + "chr1\t1000\t.\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL;END=2000;SVLEN=-1000\tGT\t"
        + "\t".join("0/1" for _ in samples) + "\n")
    gz = tmp / f"{name}.vcf.gz"
    gz.write_bytes(subprocess.run(["bgzip", "-c", str(plain)],
                                  capture_output=True).stdout)
    subprocess.run(["tabix", "-p", "vcf", str(gz)], capture_output=True)
    return gz


def run_pipeline(tmp, env, resume):
    """Returns (completed, cached), or None if the run failed."""
    trace = tmp / f"trace_{'second' if resume else 'first'}.txt"
    cmd = ["nextflow", "-q", "run", str(MAIN), "-with-docker",
           "-with-trace", str(trace),
           "--vcfs", str(tmp / "vcfs.csv"), "--ped", str(tmp / "cohort.ped"),
           "--outdir", str(tmp / "results")]
    if resume:
        cmd.append("-resume")
    # cwd=tmp so the cache DB and work dir live with the fixture, not in the repo.
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp), env=env)
    out = r.stdout + r.stderr
    which = "second" if resume else "first"
    stats = trace_stats(trace)
    if stats is None:
        check(f"the {which} run reports a task tally", False,
              out.strip().splitlines()[-1] if out.strip() else "no output")
        return None
    completed, failed, cached = stats
    if failed:
        errors = [l for l in out.splitlines() if l.startswith("[ERROR]")]
        check(f"the {which} run succeeds", False, "; ".join(errors[:2]))
        return None
    return completed, cached


def main():
    for tool in ("nextflow", "bgzip", "tabix", "docker"):
        if not shutil.which(tool):
            print(f"SKIP: {tool} not available")
            return 0

    env = dict(os.environ, NXF_SYNTAX_PARSER="v2", NXF_ANSI_LOG="false")

    # ONE directory for both runs. See the module docstring: a fresh one per run changes
    # every input file's mtime and forces a total cache miss regardless of the bug.
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # TWO sample sets, TWO callers each, so axis A merges twice and promote_axis_a runs
        # more than once. That multiplicity is the point: promote_axis_a used to take its
        # label from one channel and its VCF from a process output, which emits in task
        # COMPLETION order -- so the pairing, and therefore the task hash, varied run to run
        # and it was the one step -resume could never match. With a single group there is
        # only one item and nothing to permute, so the old single-input fixture could not
        # have caught it.
        groups = (("GA", ["S1", "S2"]), ("GB", ["S3", "S4"]))
        rows = []
        for group, samples in groups:
            for caller in ("manta", "smoove"):
                gz = write_vcf(tmp, f"{group}.{caller}", samples)
                rows.append(f"{group},{caller},false,{gz},{gz}.tbi\n")
        (tmp / "vcfs.csv").write_text("sample_set,caller,joint,vcf,tbi\n" + "".join(rows))
        (tmp / "cohort.ped").write_text(
            "".join(f"FAM_{s}\t{s}\t0\t0\t1\t0\n"
                    for _, samples in groups for s in samples))

        first = run_pipeline(tmp, env, resume=False)
        if first is None:
            return 1
        completed_first, _ = first
        check("the first run does some work", completed_first > 0,
              f"completed={completed_first}")

        second = run_pipeline(tmp, env, resume=True)
        if second is None:
            return 1
        completed_second, cached_second = second

        # The real assertion. `cached > 0` alone would pass on a pipeline that cached one
        # task out of nine, which is the shape a partial regression takes.
        check("the second run caches every task", cached_second == completed_first,
              f"cached={cached_second} of {completed_first}")
        check("the second run re-runs nothing", completed_second == 0,
              f"completed={completed_second}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: " + ", ".join(failures))
        return 1
    print("-resume caches the whole graph on an unchanged rerun")
    return 0


if __name__ == "__main__":
    sys.exit(main())
