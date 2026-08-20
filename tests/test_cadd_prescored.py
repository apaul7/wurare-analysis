#!/usr/bin/env python3
"""Regression test for the prescored CADD top-up (annotate_snps).

    python3 tests/test_cadd_prescored.py

Supplying --precomputed_cadd WITHOUT --skip_spliceai_squirls must prescore, not skip:
variants found in the table take its scores, only the remainder goes to CADD.sh, and the
merged table carries both. The passthrough pairing (table + skip flag) is covered by
test_skip_params.py; what is under test here is the top-up wiring in cadd.nf's prescored
branch, which lives in Nextflow channels and cannot be exercised any other way.

Two layers, because the full pipeline cannot COMPLETE a non-skip run locally (see
test_skip_params' module docstring: run_squirls hardcodes /usr/bin/java), and when the
splice branch dies its empty collections abort the whole run, racing the CADD branch:

  harness   a minimal workflow written here that includes `cadd` from the real module and
            feeds it real per-interval VCFs -- deterministic end to end. This is where the
            data-flow claims are asserted: prescored rows copied, uncovered variants
            CADD-scored, the fully-covered interval surviving run_cadd's empty-input guard
            (the fake CADD.sh fails on an empty VCF, like the real one), and the merged
            table single-headered, position-sorted, chr-prefixed and tabix-able.
  pipeline  the real main.nf, asserted only for wiring: a table without the skip flag must
            route into extract_prescored (submission is printed before the splice branch
            can die), and must still demand --cadd_data_dir at startup.

Needs nextflow, and a python3 on PATH with pysam (the real add_cadd_scores.py runs inside
the task). Skips cleanly without either.
"""

import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from test_skip_params import FAKES, build_env, check, failures, run_pipeline, write_exec

try:
    import pysam
except ImportError:
    print("SKIP: pysam not installed")
    sys.exit(0)

ROOT = Path(__file__).parent.parent
CADD_NF = ROOT / "nf" / "annotate_snps" / "modules" / "annotations" / "cadd.nf"

CADD_HEADER = ("##CADD GRCh38-v1.6 (c) University of Washington, Hudson-Alpha Institute "
               "for Biotechnology and Berlin Institute of Health 2013-2020. "
               "All rights reserved.\n"
               "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n")

VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

# The scores tell the rows apart in the merged table: prescored fixture rows carry
# 0.700000/7.700, the fake CADD.sh stamps 0.5/5.5 on everything it scores.
PRESCORED = "0.700000\t7.700"
CADD_SCORED = "0.5\t5.5"

# Fails on a variant-free input, like the real CADD.sh -- that is what proves run_cadd's
# empty-input guard fires for an interval the table fully covers.
FAKE_CADD = r"""
prev=""; out=""; in=""
for a in "$@"; do [ "$prev" = "-o" ] && out="$a"; prev="$a"; in="$a"; done
if ! grep -qv '^#' "$in"; then echo "fake CADD.sh: empty VCF" >&2; exit 1; fi
{ printf '##CADD GRCh38-v1.6 fake\n#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n'
  grep -v '^#' "$in" | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$4,$5,"0.5","5.5"}'
} | gzip -c > "$out"
exit 0
"""

# The smallest workflow that drives the real cadd module: two fixed intervals, a prescored
# pair, no publish -- the merged table is read out of merge_cadd's work dir by filename,
# which nothing else in this workflow produces.
HARNESS_NF = f"""nextflow.enable.types = true
include {{ cadd }} from '{CADD_NF}'

workflow {{
    main:
    vcfs = channel.of(
        record(interval: 'chr1', tag: 'chr1',
               vcf: file(params.vcf1), tbi: file(params.tbi1)),
        record(interval: 'chr2', tag: 'chr2',
               vcf: file(params.vcf2), tbi: file(params.tbi2))
    )
    cadd(vcfs, file(params.cadd_data_dir),
         [file(params.prescored), file(params.prescored_tbi)])
}}
"""


def write_vcf(path, rows):
    plain = path.with_suffix("")
    plain.write_text(VCF_HEADER + "".join(r + "\n" for r in rows))
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def write_table(path, rows):
    plain = path.with_suffix("")
    plain.write_text(CADD_HEADER + "".join(r + "\n" for r in rows))
    pysam.tabix_compress(str(plain), str(path), force=True)
    # tabix -s 1 -b 2 -e 2, in pysam's 0-based column indices.
    pysam.tabix_index(str(path), seq_col=0, start_col=1, end_col=1, force=True)


def run_harness(tmp, binf, table, tag):
    """Runs the harness workflow in its own directory; returns (trace_names, out, rows).

    `rows` is the merged table's lines, or None if merge_cadd never produced one.
    """
    import os
    rundir = tmp / f"harness_{tag}"
    rundir.mkdir()
    (rundir / "test.config").write_text(
        "process {\n    resourceLimits = [ memory: 2.GB, cpus: 2 ]\n}\n")
    env = dict(os.environ, NXF_SYNTAX_PARSER="v2",
               PATH=f"{binf}:{os.environ['PATH']}")
    trace = rundir / "trace.txt"
    r = subprocess.run(
        ["nextflow", "run", str(tmp / "cadd_harness.nf"),
         "-c", "test.config", "-with-trace", str(trace),
         "--vcf1", str(tmp / "chr1.vcf.gz"), "--tbi1", str(tmp / "chr1.vcf.gz.tbi"),
         "--vcf2", str(tmp / "chr2.vcf.gz"), "--tbi2", str(tmp / "chr2.vcf.gz.tbi"),
         "--cadd_data_dir", str(tmp / "cadd_dir"),
         "--prescored", str(tmp / table), "--prescored_tbi", str(tmp / f"{table}.tbi")],
        capture_output=True, text=True, cwd=str(rundir), env=env)
    out = r.stdout + r.stderr

    names = set()
    if trace.is_file():
        for line in trace.read_text().splitlines()[1:]:
            cols = line.split("\t")
            if len(cols) > 3:
                names.add(cols[3].split(":")[-1].split(" (")[0])

    rows = None
    for p in sorted((rundir / "work").rglob("caddv1.6.out.tsv.gz")):
        if p.name == "caddv1.6.out.tsv.gz" and p.is_file():
            with gzip.open(p, "rt") as fh:
                rows = fh.read().splitlines()
            break
    return names, out, rows


def assert_merged(rows, label):
    body = [r for r in rows if not r.startswith("#")]
    check(f"[{label}] exactly one column-header line",
          sum(1 for r in rows if r.startswith("#Chrom")) == 1, str(rows[:4]))
    check(f"[{label}] a ##CADD version line survives the merge",
          any(r.startswith("##CADD GRCh38-v1.6") for r in rows), str(rows[:2]))
    check(f"[{label}] the covered chr1 variant keeps the table's score",
          f"chr1\t100\tA\tG\t{PRESCORED}" in body, str(body))
    check(f"[{label}] the uncovered chr1 variant was scored by CADD.sh",
          f"chr1\t200\tC\tT\t{CADD_SCORED}" in body, str(body))
    check(f"[{label}] the fully-covered chr2 interval keeps its table row "
          "(and its run_cadd survived an empty VCF)",
          f"chr2\t300\tG\tA\t{PRESCORED}" in body, str(body))
    check(f"[{label}] no variant is scored twice", len(body) == 3, str(body))
    keys = [(r.split("\t")[0], int(r.split("\t")[1])) for r in body]
    check(f"[{label}] the body is position-sorted, so the real tabix would accept it",
          keys == sorted(keys), str(keys))


def main():
    if not shutil.which("nextflow"):
        print("SKIP: nextflow not on PATH")
        return 0
    probe = subprocess.run(  # the task-side python must see pysam, not just this one
        ["python3", "-c", "import pysam"], capture_output=True)
    if probe.returncode != 0:
        print("SKIP: python3 on PATH has no pysam (add_cadd_scores.py runs inside a task)")
        return 0

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # First: it creates cadd_dir/ and the fake tools the full-pipeline runs need; the
        # harness fixtures below overwrite CADD.sh inside it.
        binf_full = build_env(tmp)

        # --- harness fixtures -----------------------------------------------------------
        # chr1 holds a covered and an uncovered variant; chr2 is covered entirely, which is
        # the empty-input-guard case. Real bgzipped VCFs: extract_prescored reads and
        # rewrites them with pysam.
        (tmp / "cadd_harness.nf").write_text(HARNESS_NF)
        # Nextflow resolves bin/ against the ENTRY script's project dir -- the harness's,
        # not the module's -- so the real script is staged next to it.
        (tmp / "bin").mkdir()
        shutil.copy(ROOT / "nf" / "annotate_snps" / "bin" / "add_cadd_scores.py",
                    tmp / "bin" / "add_cadd_scores.py")
        write_vcf(tmp / "chr1.vcf.gz", ["chr1\t100\t.\tA\tG\t.\t.\t.",
                                        "chr1\t200\t.\tC\tT\t.\t.\t."])
        write_vcf(tmp / "chr2.vcf.gz", ["chr2\t300\t.\tG\tA\t.\t.\t."])
        write_table(tmp / "prescored_chr.tsv.gz",
                    [f"chr1\t100\tA\tG\t{PRESCORED}", f"chr2\t300\tG\tA\t{PRESCORED}"])
        write_table(tmp / "prescored_native.tsv.gz",
                    [f"1\t100\tA\tG\t{PRESCORED}", f"2\t300\tG\tA\t{PRESCORED}"])

        binf = tmp / "fakebin_harness"
        binf.mkdir()
        for name in ("bgzip", "tabix"):
            write_exec(binf / name, FAKES[name])
        write_exec(tmp / "cadd_dir" / "CADD.sh", FAKE_CADD)

        # --- top-up through the real cadd module ----------------------------------------
        print("a callset-named table prescores; CADD runs only for what it does not cover")
        names, out, rows = run_harness(tmp, binf, "prescored_chr.tsv.gz", "chr")
        for p in ("extract_prescored", "run_cadd", "merge_cadd"):
            check(f"{p} completes", p in names,
                  f"trace={sorted(names)}\n...{out.strip()[-500:]}")
        check("the merged table exists", rows is not None, f"...{out.strip()[-500:]}")
        if rows is not None:
            assert_merged(rows, "chr")

        print("a native-named table (raw CADD.sh output) prescores the same callset")
        _names, out, rows = run_harness(tmp, binf, "prescored_native.tsv.gz", "native")
        check("the merged table exists", rows is not None, f"...{out.strip()[-500:]}")
        if rows is not None:
            assert_merged(rows, "native")

        # --- the real main.nf routes a no-skip table into the top-up --------------------
        # The full non-skip run cannot complete locally (module docstring), so this reads
        # the submitted-process names off stdout, which appear at submission time --
        # extract_prescored is submitted the moment split_vcf completes, before the splice
        # branch's failure can tear the run down.
        print("main.nf: a table without the skip flag goes to extract_prescored")
        write_table(tmp / "topup.tsv.gz", [f"chr1\t100\tA\tG\t{PRESCORED}"])
        _trace, named, out = run_pipeline(tmp, binf_full, [
            "--precomputed_cadd", "topup.tsv.gz",
            "--precomputed_cadd_tbi", "topup.tsv.gz.tbi"], "topup")
        check("extract_prescored is submitted", "extract_prescored" in named,
              f"named={sorted(named)}\n...{out.strip()[-400:]}")

        print("main.nf: a table without the skip flag still requires cadd_data_dir")
        trace, _named, out = run_pipeline(tmp, binf_full, [
            "--precomputed_cadd", "topup.tsv.gz",
            "--precomputed_cadd_tbi", "topup.tsv.gz.tbi"], "nodatadir",
            skippable_resources=False)
        check("the startup error names cadd_data_dir",
              "cadd_data_dir" in out and "required" in out, f"...{out.strip()[-300:]}")
        check("and nothing was submitted", not trace, f"trace={sorted(trace)}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("prescored CADD top-up scores only what the table does not cover")
    return 0


if __name__ == "__main__":
    sys.exit(main())
