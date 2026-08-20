#!/usr/bin/env python3
"""End-to-end regression test for `-resume` in annotate_snps.

    python tests/test_resume_snps.py

The sibling test (test_resume.py) guards the same defect in annotate_svs, where it was found
and fixed first. annotate_snps had the identical bug and kept it for months longer: twelve
processes took a Nextflow `record(...)` in their `input:` block.

Nextflow has no value-based hash for a `record`, so its task hasher falls back to
Object.toString(), which embeds the JVM identity hash:

    9cb09c67... [Main.Scalars] Main.Scalars@217bf527

That is a fresh allocation address every run, so the task hash changed on every run and could
never match a cache entry. Records are fine in `output:` blocks and fine in channel plumbing
-- only `input:` is hashed. The processes therefore take interval/tag/vcf/tbi scalars and the
workflows destructure at the call site with multiMap.

Why this matters more here than in annotate_svs: CADD, SpliceAI and SQUIRLS are the expensive
half of a WGS run. A run that died in table_annovar re-scored the entire genome.

WHY THE FIXTURE IS BUILT ONCE AND BOTH RUNS SHARE IT. Nextflow's default cache mode hashes
each input file by path, size and LAST-MODIFIED. Building the fixture per invocation gives
new paths and new mtimes, so `cached=0` regardless of whether the bug is present, and the
test would fail forever for a reason unrelated to what it checks. Both runs also launch from
inside the fixture directory so `.nextflow/` (which holds the cache DB) and `work/` are
cleaned up with it rather than left in the working tree.

Needs `nextflow`. Skips cleanly without it.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from test_skip_params import ROOT, build_env

MAIN = ROOT / "nf" / "annotate_snps" / "main.nf"

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run_pipeline(tmp, binf, resume):
    """Returns (completed, cached), or None if the run failed."""
    env = dict(os.environ, NXF_SYNTAX_PARSER="v2",
               PATH=f"{binf}:{os.environ['PATH']}")
    cmd = ["nextflow", "run", str(MAIN),
           "-c", "test.config",
           "--vcf", "in.vcf.gz", "--tbi", "in.vcf.gz.tbi",
           "--cohort", "TEST", "--data_type", "wgs",
           "--reference.fa", "ref.fa", "--reference.fai", "ref.fa.fai",
           "--reference.dict", "ref.dict",
           "--annovar_dir", "annovar_dir",
           "--annovar_splice_scores_script", "scores.py",
           "--omim_xref", "gene_xref.txt",
           # The skip path, so the run is quick and deterministic. It still crosses
           # split_vcf, the whole ANNOVAR chain and build_cadd_humandb -- eight of the
           # twelve processes that carried the defect.
           "--skip_spliceai_squirls", "true",
           "--precomputed_cadd", "cadd.tsv.gz",
           "--precomputed_cadd_tbi", "cadd.tsv.gz.tbi",
           "--outdir", "results"]
    if resume:
        cmd.append("-resume")
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp), env=env)
    out = r.stdout + r.stderr
    m = re.search(r"completed=(\d+) failed=(\d+) cached=(\d+)", out)
    if not m:
        print(f"  FAIL  no run summary (resume={resume})\n{out.strip()[-800:]}")
        failures.append("run summary")
        return None
    completed, failed, cached = (int(m.group(i)) for i in (1, 2, 3))
    if failed:
        print(f"  FAIL  {failed} task(s) failed (resume={resume})\n{out.strip()[-800:]}")
        failures.append("failed tasks")
        return None
    return completed, cached


def main():
    if not shutil.which("nextflow"):
        print("SKIP: nextflow not on PATH")
        return 0

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        binf = build_env(tmp)

        first = run_pipeline(tmp, binf, resume=False)
        if first is None:
            return 1
        completed_first, _ = first
        check("the first run does some work", completed_first > 0,
              f"completed={completed_first}")

        second = run_pipeline(tmp, binf, resume=True)
        if second is None:
            return 1
        completed_second, cached_second = second

        # `cached > 0` alone would pass on a pipeline that cached one task out of twelve,
        # which is exactly the shape a partial regression takes -- and the shape this bug
        # would take if some but not all record inputs were destructured.
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
