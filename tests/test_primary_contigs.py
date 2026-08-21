#!/usr/bin/env python3
"""annotate_snps fans out over primary contigs only.

    python tests/test_primary_contigs.py

A real run (DRAGEN hg38) died in split_vcf on `bcftools view -r "HLA-DQB1*02:01:01"`:
region strings are parsed by splitting at colons, so a contig whose NAME contains colons
is misread as contig + coordinates. The same raw interval string later feeds pysam
fetch(region=...) in add_scores.py / add_cadd_scores.py, so quoting can never fix it.
The pipeline instead filters the interval list to primary contigs (chr1-22/X/Y/M/MT,
`chr` prefix optional) in subworkflows/annotate_snps.nf, warning per skipped contig --
defensible because CADD/SpliceAI/ANNOVAR carry no data for alt/decoy/HLA contigs anyway.

This runs the real DAG with the fake toolchain from test_skip_params, but with a fake
`bcftools query` that reports a mixed contig list, and asserts the fan-out happened for
exactly the primary names: split outputs exist for chr1/chrM/MT, none exist for the HLA
or alt contig, and each drop was warned about.

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

from test_skip_params import FAKES, ROOT, build_env, trace_stats, write_exec

MAIN = ROOT / "nf" / "annotate_snps" / "main.nf"

KEPT = ["chr1", "chrM", "MT"]
DROPPED = ["HLA-DQB1*02:01:01", "chr1_KI270766v1_alt"]

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

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        binf = build_env(tmp)

        # Same fake bcftools, but `query` (what list_chroms calls) reports a DRAGEN-shaped
        # contig mix instead of chr1/chr2. The replace target is pinned to the fixture's
        # current text: if test_skip_params rewords it, this fails loudly here rather than
        # silently testing the wrong contig list.
        contigs = "\\n".join(KEPT + DROPPED) + "\\n"
        body = FAKES["bcftools"].replace(
            r"printf 'chr1\nchr2\n'", f"printf '{contigs}'")
        if body == FAKES["bcftools"]:
            print("  FAIL  could not patch the fake bcftools query branch")
            return 1
        write_exec(binf / "bcftools", body)

        env = dict(os.environ, NXF_SYNTAX_PARSER="v2", NXF_ANSI_LOG="false",
                   PATH=f"{binf}:{os.environ['PATH']}")
        # The skip path, as in test_resume_snps: quick, and it still crosses list_chroms,
        # the interval filter under test, split_vcf and the whole ANNOVAR chain.
        cmd = ["nextflow", "run", str(MAIN),
               "-c", "test.config",
               "-with-trace", "trace.txt",
               "--vcf", "in.vcf.gz", "--tbi", "in.vcf.gz.tbi",
               "--cohort", "TEST", "--data_type", "wgs",
               "--reference.fa", "ref.fa", "--reference.fai", "ref.fa.fai",
               "--reference.dict", "ref.dict",
               "--annovar_dir", "annovar_dir",
               "--annovar_splice_scores_script", "scores.py",
               "--omim_xref", "gene_xref.txt",
               "--skip_spliceai_squirls", "true",
               "--precomputed_cadd", "cadd.tsv.gz",
               "--precomputed_cadd_tbi", "cadd.tsv.gz.tbi",
               "--outdir", "results"]
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp), env=env)
        out = r.stdout + r.stderr

        stats = trace_stats(tmp / "trace.txt")
        if stats is None:
            print(f"  FAIL  no trace written\n{out.strip()[-800:]}")
            return 1
        _completed, failed, _cached = stats
        check("the run completes with no failed tasks", failed == 0,
              f"failed={failed}\n{out.strip()[-800:]}")

        work = tmp / "work"
        for c in KEPT:
            check(f"split output exists for {c}",
                  any(work.rglob(f"{c}.out.vcf.gz")))
        for c in DROPPED:
            tag = re.sub(r"[^A-Za-z0-9._-]", "_", c)
            check(f"no split output for {c}",
                  not any(work.rglob(f"{tag}.out.vcf.gz")))
            check(f"the drop of {c} is warned about",
                  f"skipping non-primary contig '{c}'" in out)

    print()
    if failures:
        print(f"{len(failures)} FAILED: " + ", ".join(failures))
        return 1
    print("annotate_snps fans out over primary contigs only, warning per dropped contig")
    return 0


if __name__ == "__main__":
    sys.exit(main())
