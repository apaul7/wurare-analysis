#!/usr/bin/env python3
"""Regression test for annotate_snps' re-annotation skip params.

    python3 tests/test_skip_params.py

Like test_grouping_e2e.py this runs the real pipeline, because what it guards lives in
Nextflow channel wiring and cannot be reproduced any other way. Unlike that test it needs no
containers: it puts fake tools on PATH -- bcftools, tabix, bedtools, CADD.sh and the four
ANNOVAR perl scripts -- so the skip path runs end to end on a laptop in seconds.

What is asserted, and why each direction is needed:

  skips ON   CADD/SpliceAI/SQUIRLS must NOT run, and the ANNOVAR chain MUST still complete.
             A skip that also severed ANNOVAR's input would satisfy "did it skip?" while
             producing nothing at all.
  skips OFF  the same stages MUST run -- otherwise the absence assertions above pass
             vacuously against a pipeline that is simply broken for unrelated reasons. The
             first version of this test had no control and reported a clean pass while the
             run was dying before it ever reached CADD.

Two different observables, deliberately:

  skip mode  runs to completion, so the -with-trace TSV is complete and is the source of
             truth. Trace records only tasks that reached an exit status.
  no-skip    cannot complete locally: run_squirls hardcodes /usr/bin/java, an absolute path
             no PATH shim can intercept, and add_precomputed needs pysam. So the control
             reads the qualified process names Nextflow prints to stdout, which appear when
             a task is submitted regardless of whether it later fails.

Why not `nextflow lint`/`inspect`: measured, both miss this class entirely. `inspect`
statically collects every included process without executing the workflow body, so a
conditional that skips a stage looks identical to one that does not -- an emit-shadowing
bug passed both and was caught only by running.

Also covered: the paired-param guard (precomputed_cadd without its .tbi must fail at startup,
not inside a task), and that build_cadd_humandb still runs under precomputed_cadd -- supplying
a scored table skips the scoring, not the humandb reformat/index ANNOVAR needs either way.
"""

import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
MAIN = ROOT / "nf" / "annotate_snps" / "main.nf"

# Stages that the skip params must remove from the run. extract_prescored belongs here even
# though the no-skip control cannot reach it (it needs a --precomputed_cadd): with the table
# AND the skip flag the run is a passthrough, and extracting would prove the passthrough
# branch was not taken. Its positive direction is test_cadd_prescored.py's job.
SKIPPABLE = {"run_cadd", "merge_cadd", "extract_prescored", "run_spliceai",
             "add_precomputed", "concat", "run_squirls", "compress_squirls", "merge_vcfs"}

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


# --- fake tools ---------------------------------------------------------------------------
# Each is the smallest thing that produces the output files the process declares. They are
# not simulations: nothing here checks flags beyond what it must to place output correctly.

FAKES = {
    "bcftools": r"""#!/bin/sh
sub="$1"; shift
out=""; prev=""; last=""
for a in "$@"; do
  [ "$prev" = "-o" ] && out="$a"
  prev="$a"; last="$a"
done
hdr='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
row='chr1\t100\t.\tA\tG\t99\tPASS\t.\tGT\t0/1\n'
case "$sub" in
  query)     printf 'chr1\nchr2\n' ;;
  view|norm|concat|sort)
             [ -n "$out" ] && printf "$hdr$row" > "$out" ;;
  index)     : > "$last.tbi" ;;
  --version) echo "bcftools 1.99-fake" ;;
esac
exit 0
""",
    "tabix": r"""#!/bin/sh
case "$1" in --version) echo "tabix (htslib) 1.99-fake"; exit 0 ;; esac
for a in "$@"; do last="$a"; done
: > "$last.tbi"
exit 0
""",
    # bgzip must really compress: merge_cadd's output is read back with zgrep.
    "bgzip": r"""#!/bin/sh
if [ "$1" = "-c" ]; then shift; gzip -c "$@"; else gzip -f "$@"; fi
exit 0
""",
    "bedtools": r"""#!/bin/sh
case "$1" in --version) echo "bedtools v2.99-fake"; exit 0 ;; esac
prev=""; in=""
for a in "$@"; do [ "$prev" = "-i" ] && in="$a"; prev="$a"; done
[ -n "$in" ] && cat "$in"
exit 0
""",
}

# ANNOVAR perl scripts. Bare invocation prints a version banner and exits non-zero, which is
# what annovar_version() in annovar.nf probes for -- faking that keeps versions.yml honest.
ANNOVAR_BANNER = r"""
if [ $# -eq 0 ]; then echo "Version: \$Date: 2020-06-08 00:00:00 -0400 (Mon, 8 Jun 2020) \$"; exit 1; fi
"""

ANNOVAR_FAKES = {
    # index_annovar.pl --filetype A --outfile humandb/hg38_CADDv1.6.txt <in>
    "index_annovar.pl": ANNOVAR_BANNER + r"""
prev=""; out=""
for a in "$@"; do [ "$prev" = "--outfile" ] && out="$a"; prev="$a"; done
: > "$out"; : > "$out.idx"
exit 0
""",
    # convert2annovar.pl -format vcf4 <vcf> -outfile sample --includeinfo -allsample
    "convert2annovar.pl": ANNOVAR_BANNER + r"""
prev=""; out="sample"
for a in "$@"; do [ "$prev" = "-outfile" ] && out="$a"; prev="$a"; done
printf 'chr1\t100\t100\tA\tG\thet\t99\t20\n' > "$out.S1.avinput"
exit 0
""",
    # variants_reduction.pl <avinput> ./humandb ... --outfile <prefix>
    "variants_reduction.pl": ANNOVAR_BANNER + r"""
prev=""; out=""
for a in "$@"; do [ "$prev" = "--outfile" ] && out="$a"; prev="$a"; done
printf 'chr1\t100\t100\tA\tG\n' > "$out.step2.varlist"
exit 0
""",
    # table_annovar.pl <varlist> ./humandb/ ... --out <prefix>
    "table_annovar.pl": ANNOVAR_BANNER + r"""
prev=""; out=""
for a in "$@"; do [ "$prev" = "--out" ] && out="$a"; prev="$a"; done
printf 'Chr\tStart\tEnd\tRef\tAlt\tCADD_phred\n' > "$out.hg38_multianno.txt"
printf 'chr1\t100\t100\tA\tG\t9.566\n' >> "$out.hg38_multianno.txt"
exit 0
""",
}


def write_exec(path, body):
    path.write_text(body if body.startswith("#!") else "#!/bin/sh\n" + body)
    path.chmod(0o755)


# Derived from conf/params.config's annovar_protocols rather than hardcoded, so adding a
# protocol there cannot silently leave the fixture -- and therefore the preflight -- uncovered.
# CLINVAR is the placeholder main.nf substitutes with --clinvar_date.
PARAMS_CONFIG = ROOT / "nf" / "annotate_snps" / "conf" / "params.config"
CLINVAR_DATE = "clinvar_20260627"

# Two protocols do not follow hg38_<name>.txt, and one is built at runtime.
HUMANDB_SPECIAL = {
    "1000g2015aug_all": ["hg38_ALL.sites.2015_08.txt"],
    "refGene": ["hg38_refGene.txt", "hg38_refGeneMrna.fa"],
    "CADDv1.6": [],
}


def annovar_protocols():
    """The protocol list main.nf preflights against, read from params.config."""
    text = PARAMS_CONFIG.read_text()
    m = re.search(r'annovar_protocols\s*=\s*(.+?)\n\s*annovar_operations', text, re.S)
    if not m:
        raise RuntimeError("could not find annovar_protocols in params.config")
    joined = "".join(re.findall(r'"([^"]*)"', m.group(1)))
    return joined.replace("CLINVAR", CLINVAR_DATE).split(",")


def humandb_stubs():
    """Every file the preflight requires, plus one .idx for stage_humandb's glob."""
    out = []
    for p in annovar_protocols():
        out += HUMANDB_SPECIAL.get(p, [f"hg38_{p}.txt"])
    # No individual .idx is required, but stage_humandb globs *.txt.idx and `ln -s` fails on
    # an unmatched glob under `set -euo pipefail`.
    out.append("hg38_refGene.txt.idx")
    return out


def build_env(tmp):
    """Creates every fake tool and input file. Returns the fakebin dir."""
    binf = tmp / "fakebin"
    binf.mkdir()
    for name, body in FAKES.items():
        write_exec(binf / name, body)

    # ANNOVAR install: the perl scripts, plus a humandb/ holding one file per glob that
    # stage_humandb() expands. An unmatched glob would make `ln -s` fail under `set -e`.
    annovar = tmp / "annovar_dir"
    (annovar / "humandb").mkdir(parents=True)
    for name, body in ANNOVAR_FAKES.items():
        write_exec(annovar / name, body)
    for stub in humandb_stubs():
        (annovar / "humandb" / stub).touch()

    # CADD install, for the no-skip control.
    cadd_dir = tmp / "cadd_dir"
    cadd_dir.mkdir()
    write_exec(cadd_dir / "CADD.sh", r"""
prev=""; out=""
for a in "$@"; do [ "$prev" = "-o" ] && out="$a"; prev="$a"; done
printf '#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t100\tA\tG\t0.8\t9.566\n' | gzip -c > "$out"
exit 0
""")
    write_exec(binf / "spliceai", "exit 0\n")

    # A fake somalier, for the QC stage. It fakes the TOOL, never the contract: the real
    # column names and the real relatedness semantics are pinned against the real container
    # in tests/test_somalier_assumptions.py. What is under test here is the wiring --
    # extract -> relate -> check_somalier, and the 04_qc publish targets.
    write_exec(binf / "somalier", r"""
sub="$1"; shift
case "$sub" in
  --version|version) echo "somalier version: 0.2.19-fake"; exit 0 ;;
  extract)
    d="."
    while [ $# -gt 0 ]; do case "$1" in -d) d="$2"; shift 2 ;; *) shift ;; esac; done
    mkdir -p "$d"; : > "$d/SAMP_M.somalier"; : > "$d/SAMP_F.somalier"; exit 0 ;;
  relate)
    # Full ped-shaped samples.tsv, like the real tool writes: the first six columns are what
    # somalier_relate cuts into the draft ped. Parents are -9 -- somalier's "unknown" -- so
    # the draft-ped normalization to 0 is exercised. Sexes agree with the X/Y depths.
    # The trailing "# relate-args:" comment records what this fake was invoked with, so the
    # test can assert --infer routing; every consumer skips comment/short rows.
    o="cohort"; args="$*"
    while [ $# -gt 0 ]; do case "$1" in -o) o="$2"; shift 2 ;; *) shift ;; esac; done
    printf '#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tdepth_mean\tX_depth_mean\tY_depth_mean\n' > "$o.samples.tsv"
    printf 'FAM\tSAMP_M\t-9\t-9\t1\t-9\t30\t15\t12\n' >> "$o.samples.tsv"
    printf 'FAM\tSAMP_F\t-9\t-9\t2\t-9\t30\t30\t0.3\n' >> "$o.samples.tsv"
    printf '# relate-args: %s\n' "$args" >> "$o.samples.tsv"
    printf '#sample_a\tsample_b\trelatedness\texpected_relatedness\n' > "$o.pairs.tsv"
    printf 'SAMP_M\tSAMP_F\t0.02\t-1.0\n' >> "$o.pairs.tsv"
    : > "$o.html"; exit 0 ;;
  ancestry)
    # The DOUBLED filename is the point of this fake. somalier appends a fixed
    # ".somalier-ancestry.<ext>" to its -o prefix, whose default is itself
    # "somalier-ancestry" -- so it writes somalier-ancestry.somalier-ancestry.tsv. A fake
    # that wrote the name the process DECLARES would pass against the broken code and prove
    # nothing. The real image is pinned in tests/test_somalier_assumptions.py.
    #
    # The table carries the labelled 1kg samples as well as the query ones, matching the real
    # tool: ancestry.nim writes a row per labelled sample and per query sample, with
    # given_ancestry filled for the former and empty for the latter.
    q=""; after=0
    for a in "$@"; do
      if [ "$a" = "++" ]; then after=1; continue; fi
      [ "$after" = 1 ] && q="$q $a"
    done
    out=somalier-ancestry.somalier-ancestry.tsv
    printf '#sample_id\tpredicted_ancestry\tgiven_ancestry\n' > "$out"
    printf 'HG00096\tEUR\tEUR\n' >> "$out"
    printf 'HG00097\tAFR\tAFR\n' >> "$out"
    for f in $q; do
      s=$(basename "$f" .somalier)
      printf '%s\tEUR\t\n' "$s" >> "$out"
    done
    : > somalier-ancestry.somalier-ancestry.html
    exit 0 ;;
esac
exit 0
""")
    (tmp / "sites.vcf.gz").touch()
    # Ancestry fixtures. The 1kg directory needs real files: somalier_ancestry globs
    # <dir>/*.somalier, and under `set -euo pipefail` an unmatched glob would reach the tool
    # as a literal word.
    (tmp / "ancestry-labels.tsv").write_text("sample\tlabel\nHG00096\tEUR\nHG00097\tAFR\n")
    kg_dir = tmp / "1kg"
    kg_dir.mkdir()
    for bg in ("HG00096", "HG00097"):
        (kg_dir / f"{bg}.somalier").touch()
    # SAMP_F's PED sex is wrong on purpose: the QC stage must override it from the data.
    (tmp / "cohort.ped").write_text("FAM\tSAMP_M\t0\t0\t1\t2\nFAM\tSAMP_F\t0\t0\t1\t1\n")

    # add_splice_scores runs this with --input/--output; copying is enough.
    write_exec(tmp / "scores.py", r"""#!/usr/bin/env python3
import sys
a = dict(zip(sys.argv[1::2], sys.argv[2::2]))
open(a["--output"], "w").write(open(a["--input"]).read())
""")

    # A REAL gzip: build_cadd_humandb reads this with zgrep, so an empty file aborts the task
    # before anything under test runs. That is what made the first version of this test
    # report failures that had nothing to do with the skip params.
    with gzip.open(tmp / "cadd.tsv.gz", "wt") as fh:
        fh.write("#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n")
        fh.write("chr1\t100\tA\tG\t0.8\t9.566\n")
    (tmp / "cadd.tsv.gz.tbi").touch()

    for n in ("in.vcf.gz", "in.vcf.gz.tbi", "ref.fa", "ref.fa.fai", "ref.dict",
              "squirls.yml", "jannovar.ser", "splice.vcf.gz", "splice.vcf.gz.tbi",
              "gene_xref.txt"):
        (tmp / n).touch()

    # The pipeline's real memory requests (build_cadd_humandb asks 20 GB, reduce_variants
    # and table_annovar 50 GB) exceed any laptop, and Nextflow refuses to start a task it
    # cannot satisfy -- "Process requirement exceeds available memory", with the task never
    # running and an empty work dir. resourceLimits caps the request without touching the
    # declarations, which is what compute2.config does for the cluster.
    (tmp / "test.config").write_text(
        "process {\n    resourceLimits = [ memory: 2.GB, cpus: 2 ]\n}\n")
    return binf


def run_pipeline(tmp, binf, extra_args, tag, skippable_resources=True):
    """Returns (processes_in_trace, processes_named_on_stdout, combined_output).

    `skippable_resources=False` omits the five params that only the CADD, SQUIRLS and
    SpliceAI stages read. That is what the README's re-annotation command actually passes,
    and it used to fail at startup because all eight site-local resources were demanded
    unconditionally -- so the documented command could not be copy-pasted.
    """
    trace = tmp / f"trace_{tag}.txt"
    env = dict(os.environ, NXF_SYNTAX_PARSER="v2",
               PATH=f"{binf}:{os.environ['PATH']}")
    cmd = ["nextflow", "run", str(MAIN),
           "-c", "test.config",
           "-with-trace", str(trace),
           "--vcf", "in.vcf.gz", "--tbi", "in.vcf.gz.tbi",
           "--cohort", "TEST", "--data_type", "wgs",
           "--reference.fa", "ref.fa", "--reference.fai", "ref.fa.fai",
           "--reference.dict", "ref.dict"]
    if skippable_resources:
        cmd += ["--cadd_data_dir", "cadd_dir",
                "--squirls_config", "squirls.yml",
                "--squirls_jannovar_model", "jannovar.ser",
                "--spliceai_precomputed_scores", "splice.vcf.gz",
                "--spliceai_precomputed_tbi", "splice.vcf.gz.tbi"]
    cmd += ["--annovar_dir", "annovar_dir",
            "--annovar_splice_scores_script", "scores.py",
            "--omim_xref", "gene_xref.txt",
            "--outdir", f"results_{tag}"] + extra_args
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp), env=env)
    out = r.stdout + r.stderr

    in_trace = set()
    if trace.is_file():
        for line in trace.read_text().splitlines()[1:]:
            cols = line.split("\t")
            if len(cols) > 3:
                # "annotate_snps:cadd:run_cadd (1)" -> "run_cadd"
                in_trace.add(cols[3].split(":")[-1].split(" (")[0])
    # Same shape, but scraped from console output, which names a process when it is
    # submitted -- so it survives a run that dies partway.
    named = {m.split(":")[-1] for m in re.findall(r"annotate_snps:[A-Za-z0-9_:]+", out)}
    return in_trace, named, out


def relate_args(samples_tsv):
    """The argv the fake somalier `relate` recorded into its samples.tsv."""
    if not samples_tsv.is_file():
        return ""
    for line in samples_tsv.read_text().splitlines():
        if line.startswith("# relate-args:"):
            return line.split(":", 1)[1]
    return ""


def main():
    if not shutil.which("nextflow"):
        print("SKIP: nextflow not on PATH")
        return 0

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        binf = build_env(tmp)

        # --- control ------------------------------------------------------------------
        # Cannot complete locally (run_squirls hardcodes /usr/bin/java), so this reads the
        # submitted-process names off stdout rather than the trace.
        print("skips OFF -- the stages under test must actually run")
        _tr, named, out = run_pipeline(tmp, binf, [], "noskip")
        check("the run reaches the branch point",
              "split_vcf" in named, f"named={sorted(named)}")
        check("a skippable stage is submitted",
              bool(SKIPPABLE & named),
              f"none of {sorted(SKIPPABLE)} in {sorted(named)} -- the control is broken, so "
              "the absence assertions below would prove nothing")

        # --- both skips on --------------------------------------------------------------
        print("skips ON -- stages must not run, and ANNOVAR must still complete")
        trace, named, out = run_pipeline(tmp, binf, [
            "--skip_spliceai_squirls", "true",
            "--precomputed_cadd", "cadd.tsv.gz",
            "--precomputed_cadd_tbi", "cadd.tsv.gz.tbi"], "skip")
        check("the skip run completes with no failed task",
              "[FAILED]" not in out and "ERROR" not in out,
              f"...{out.strip()[-400:]}")
        check("no CADD/SpliceAI/SQUIRLS stage ran",
              not (SKIPPABLE & (trace | named)),
              f"ran {sorted(SKIPPABLE & (trace | named))} -- the skip did not take effect")
        check("split_vcf still runs (intervals still come from the input VCF)",
              "split_vcf" in trace, f"trace={sorted(trace)}")
        check("normalize_vcf runs, so ANNOVAR is still wired to the split VCFs",
              "normalize_vcf" in trace,
              f"trace={sorted(trace)} -- the skip severed ANNOVAR's input instead of "
              "rerouting it, which an absence-only assertion cannot catch")
        check("build_cadd_humandb still runs on the supplied table",
              "build_cadd_humandb" in trace,
              f"trace={sorted(trace)} -- supplying a scored table skips scoring, not the "
              "humandb reformat/index ANNOVAR needs regardless")
        check("the ANNOVAR chain reaches merge_annovar",
              "merge_annovar" in trace, f"trace={sorted(trace)}")

        # --- the README's command, exactly as documented -----------------------------------
        # The re-annotation example passes neither --cadd_data_dir nor the SQUIRLS/SpliceAI
        # paths, because nothing that runs reads them. All eight site-local resources used to
        # be required unconditionally, so the documented command died at startup naming five
        # params -- and each then had to EXIST on disk, not merely be set.
        print("the documented re-annotation command runs without the skipped stages' inputs")
        trace2, _named2, out2 = run_pipeline(tmp, binf, [
            "--skip_spliceai_squirls", "true",
            "--precomputed_cadd", "cadd.tsv.gz",
            "--precomputed_cadd_tbi", "cadd.tsv.gz.tbi"], "readme",
            skippable_resources=False)
        check("it is not rejected at startup",
              "is required and has no default" not in out2,
              f"...{out2.strip()[-400:]}")
        check("and it completes",
              "[FAILED]" not in out2 and "ERROR" not in out2,
              f"...{out2.strip()[-400:]}")
        check("still producing the merged multianno TSV",
              "merge_annovar" in trace2, f"trace={sorted(trace2)}")

        tsv = tmp / "results_skip" / "02_annovar" / "TEST_wgs_20260627.hg38_multianno.tsv"
        check("a merged multianno TSV is published", tsv.is_file(), str(tsv))
        if tsv.is_file():
            check("with exactly one header line",
                  sum(1 for l in tsv.read_text().splitlines() if l.startswith("Chr")) == 1,
                  tsv.read_text()[:200])
        check("and no annotated VCF is republished",
              not (tmp / "results_skip" / "03_annotated_vcf").exists(),
              "skipping should not re-emit a copy of the VCF the caller supplied")

        # --- paired-param guard ---------------------------------------------------------
        print("precomputed_cadd without its index -- must fail at startup")
        trace, _named, out = run_pipeline(tmp, binf, ["--precomputed_cadd", "cadd.tsv.gz"],
                                          "unpaired")
        check("the error names both params",
              "precomputed_cadd_tbi" in out and "together" in out,
              f"...{out.strip()[-300:]}")
        check("and nothing was submitted",
              not trace,
              f"trace={sorted(trace)} -- the guard fired inside a task rather than at startup")

        # --- ANNOVAR humandb preflight --------------------------------------------------
        # The point of the check is that table_annovar.pl otherwise dies AFTER CADD and
        # SpliceAI. A preflight that never fires is indistinguishable from no preflight, so
        # assert on the name of the missing database, not merely on a non-zero exit.
        print("humandb preflight -- a missing database must fail at startup, by name")
        victim = tmp / "annovar_dir" / "humandb" / "hg38_varity.txt"
        victim.rename(victim.with_suffix(".txt.hidden"))
        trace, _named, out = run_pipeline(tmp, binf, [], "nodb")
        victim.with_suffix(".txt.hidden").rename(victim)
        check("the error names the missing protocol",
              "varity" in out,
              f"...{out.strip()[-400:]}")
        check("and says it was a preflight failure",
              "ANNOVAR preflight failed" in out,
              f"...{out.strip()[-400:]}")
        check("and nothing was submitted",
              not trace,
              f"trace={sorted(trace)} -- the check ran inside a task, which is the failure "
              "mode it exists to prevent")

        # refGene needs the transcript FASTA as well as the .txt -- a per-protocol special
        # case that a naive hg38_<name>.txt loop would miss.
        print("humandb preflight -- refGene's FASTA is required too")
        fa = tmp / "annovar_dir" / "humandb" / "hg38_refGeneMrna.fa"
        fa.rename(fa.with_suffix(".fa.hidden"))
        _trace, _named, out = run_pipeline(tmp, binf, [], "nofa")
        fa.with_suffix(".fa.hidden").rename(fa)
        check("the error names refGeneMrna.fa",
              "refGeneMrna.fa" in out, f"...{out.strip()[-400:]}")

        # --clinvar_date names a database, so overriding it must move what is checked. This
        # is the re-annotation path the README documents.
        print("humandb preflight -- --clinvar_date changes which database is required")
        _trace, _named, out = run_pipeline(
            tmp, binf, ["--clinvar_date", "clinvar_20991231"], "clinvarshift")
        check("the error names the overridden clinvar database",
              "clinvar_20991231" in out, f"...{out.strip()[-400:]}")

        # protocol/operation lists are paired by position; a mismatch silently applies the
        # wrong operation to every database after the gap.
        print("protocol/operation count mismatch -- must fail at startup")
        _trace, _named, out = run_pipeline(
            tmp, binf, ["--annovar_operations", "gx,r,f"], "opmismatch")
        check("the error names both counts",
              "annovar_operations" in out and "position" in out,
              f"...{out.strip()[-400:]}")

        # --- sample QC (somalier) -----------------------------------------------------
        # Optional and off by default, so the whole chain has to be shown to run when it is
        # asked for -- and to stay absent when it is not. Neither state errors on its own,
        # which is exactly why both are asserted.
        print("sample QC does not run unless it is asked for")
        trace, _named, _out = run_pipeline(tmp, binf, [], "noqc")
        check("no somalier process is submitted",
              not {p for p in trace if p.startswith("somalier")
                   or p == "check_somalier"},
              str(sorted(trace)))

        print("sample QC runs, and the ploidy table is derived from what somalier reported")
        # --skip_spliceai_squirls, because this is the one case here that has to RUN TO
        # COMPLETION rather than just reach submission: nothing is published by a failed run,
        # and add_precomputed needs pysam, which the fake environment has no way to provide.
        # Skipping it costs nothing -- the QC stage is independent of every annotation stage.
        trace, _named, out = run_pipeline(
            tmp, binf,
            ["--skip_spliceai_squirls",
             "--somalier_sites", "sites.vcf.gz", "--ped", "cohort.ped"], "qc")
        for proc in ("somalier_extract_vcf", "somalier_relate", "check_somalier"):
            check(f"{proc} ran", proc in trace, str(sorted(trace)))

        qc_dir = tmp / "results_qc" / "04_qc"
        ploidy = qc_dir / "ploidy.tsv"
        check("04_qc/ploidy.tsv is published", ploidy.is_file(),
              str(sorted(p.name for p in qc_dir.glob("*")) if qc_dir.is_dir()
                  else "no 04_qc"))
        if ploidy.is_file():
            rows = {r[0]: r[1:] for r in
                    (l.split("\t") for l in ploidy.read_text().splitlines()
                     if not l.startswith("#") and l.strip())}
            check("the hemizygous sample is called 1 X",
                  rows.get("SAMP_M", [])[:2] == ["1", "1"], str(rows))
            check("the diploid sample is called 2 X",
                  rows.get("SAMP_F", [])[:2] == ["2", "0"], str(rows))
            check("the PED's wrong sex is reported rather than believed",
                  rows.get("SAMP_F", [])[3] == "DISAGREES", str(rows))
        for name in ("cohort.inferred.ped", "TEST.samples.tsv", "TEST.pairs.tsv"):
            check(f"04_qc/{name} is published", (qc_dir / name).is_file(),
                  str(sorted(p.name for p in qc_dir.glob("*")) if qc_dir.is_dir()
                      else "no 04_qc"))
        args = relate_args(qc_dir / "TEST.samples.tsv")
        check("relate got the operator's --ped and NOT --infer -- inference must never "
              "compete with a stated pedigree",
              "--ped" in args and "--infer" not in args, repr(args))

        # The PED is optional. Without one, somalier runs with --infer and its own inferred
        # pedigree (the draft ped somalier_relate cuts from samples.tsv) stands in for the
        # operator's -- so check_somalier still runs and 04_qc still gains ploidy.tsv and a
        # cohort.inferred.ped the SV pipeline's --ped can consume. This used to be an error
        # at startup; requiring a pedigree in order to look at the data gated the check on
        # paperwork.
        print("sample QC runs WITHOUT a PED -- somalier still reports observed relatedness")
        trace, _named, out = run_pipeline(
            tmp, binf,
            ["--skip_spliceai_squirls", "--somalier_sites", "sites.vcf.gz"], "qcnoped")
        check("it is not rejected at startup",
              "params.ped is required" not in out, f"...{out.strip()[-400:]}")
        check("and the run completes",
              "[FAILED]" not in out and "ERROR" not in out, f"...{out.strip()[-400:]}")
        for proc in ("somalier_extract_vcf", "somalier_relate"):
            check(f"{proc} still ran", proc in trace, str(sorted(trace)))
        check("check_somalier still runs -- the draft ped stands in for the operator's",
              "check_somalier" in trace, str(sorted(trace)))

        # Ancestry publishes TWO tables. The full one carries the 1kg reference samples the
        # prediction was made against; the cohort one is the same table cut down to this run.
        #
        # Nothing covered this stage before, which is how it shipped declaring output files
        # somalier does not write: the tool appends a fixed ".somalier-ancestry.<ext>" to a
        # prefix that already defaults to "somalier-ancestry", so the real names are doubled
        # and every ancestry run died on "Missing output file(s)".
        print("ancestry publishes a full table and a cohort-only one")
        trace, _named, out = run_pipeline(
            tmp, binf,
            ["--skip_spliceai_squirls",
             "--somalier_sites", "sites.vcf.gz", "--ped", "cohort.ped",
             "--somalier_labels", "ancestry-labels.tsv",
             "--somalier_1kg_dir", "1kg"], "ancestry")
        check("somalier_ancestry ran", "somalier_ancestry" in trace, str(sorted(trace)))
        check("the run completes",
              "[FAILED]" not in out and "ERROR" not in out, f"...{out.strip()[-500:]}")

        anc_dir = tmp / "results_ancestry" / "04_qc"
        full = anc_dir / "somalier-ancestry.tsv"
        cohort = anc_dir / "TEST.somalier-ancestry.tsv"
        listing = (str(sorted(p.name for p in anc_dir.glob("*"))) if anc_dir.is_dir()
                   else "no 04_qc")
        check("04_qc/somalier-ancestry.tsv is published", full.is_file(), listing)
        check("04_qc/TEST.somalier-ancestry.tsv is published", cohort.is_file(), listing)
        check("the doubled name is not published", not (
            anc_dir / "somalier-ancestry.somalier-ancestry.tsv").is_file(), listing)
        check("04_qc/somalier-ancestry.html is published",
              (anc_dir / "somalier-ancestry.html").is_file(), listing)

        if full.is_file() and cohort.is_file():
            full_ids = {l.split("\t")[0] for l in full.read_text().splitlines()
                        if l.strip() and not l.startswith("#")}
            cohort_ids = {l.split("\t")[0] for l in cohort.read_text().splitlines()
                          if l.strip() and not l.startswith("#")}
            check("the full table carries the 1kg reference samples",
                  {"HG00096", "HG00097"} <= full_ids, str(sorted(full_ids)))
            check("the full table also carries this run's samples",
                  {"SAMP_M", "SAMP_F"} <= full_ids, str(sorted(full_ids)))
            check("the cohort table is exactly this run's samples",
                  cohort_ids == {"SAMP_M", "SAMP_F"}, str(sorted(cohort_ids)))
            check("the cohort table keeps its header",
                  cohort.read_text().startswith("#sample_id"),
                  cohort.read_text()[:80])

        noped_dir = tmp / "results_qcnoped" / "04_qc"
        for name in ("TEST.samples.tsv", "TEST.pairs.tsv",
                     "ploidy.tsv", "cohort.inferred.ped"):
            check(f"04_qc/{name} is still published -- the draft ped stands in for "
                  "a missing operator PED that used to gate publication",
                  (noped_dir / name).is_file(),
                  str(sorted(p.name for p in noped_dir.glob("*")) if noped_dir.is_dir()
                      else "no 04_qc"))

        args = relate_args(noped_dir / "TEST.samples.tsv")
        check("relate got --infer and no --ped -- somalier's inference is the sole source "
              "of family structure when the operator states none",
              "--infer" in args and "--ped" not in args, repr(args))

        # The published inferred ped is what an SV run's --ped would consume, so it must
        # survive prepare_svs' preflight (nf/annotate_svs/subworkflows/prepare_svs.nf):
        # >= 6 whitespace-separated columns per data line, and the sample column naming
        # this cohort. The draft ped is headerless and write_inferred_ped only reproduces
        # input lines, so no comment/header rows should appear at all.
        inferred = noped_dir / "cohort.inferred.ped"
        if inferred.is_file():
            lines = [l for l in inferred.read_text().splitlines() if l.strip()]
            check("the inferred ped carries no header or comment lines",
                  not any(l.lstrip().startswith("#") for l in lines), str(lines))
            rows = [re.split(r"\s+", l.strip()) for l in lines
                    if not l.lstrip().startswith("#")]
            check("every row has the >= 6 columns prepare_svs' preflight demands",
                  rows and all(len(r) >= 6 for r in rows), str(rows))
            check("its sample column is exactly this cohort",
                  {r[1] for r in rows if len(r) >= 2} == {"SAMP_M", "SAMP_F"}, str(rows))
            by = {r[1]: r for r in rows if len(r) >= 6}
            check("the depth cross-check wrote sex 1 for the hemizygous sample",
                  by.get("SAMP_M", [""] * 6)[4] == "1", str(by))
            check("and sex 2 for the diploid sample",
                  by.get("SAMP_F", [""] * 6)[4] == "2", str(by))
            check("somalier's -9 unknown parents were normalized to PED's 0",
                  all(r[2] == "0" and r[3] == "0" for r in rows), str(rows))

        noped_ploidy = noped_dir / "ploidy.tsv"
        if noped_ploidy.is_file():
            prows = {r[0]: r[1:] for r in
                     (l.split("\t") for l in noped_ploidy.read_text().splitlines()
                      if not l.startswith("#") and l.strip())}
            check("both samples AGREE with the draft ped -- the same depths wrote both "
                  "sides of the comparison",
                  all(prows.get(s, ["", "", "", ""])[3] == "AGREES"
                      for s in ("SAMP_M", "SAMP_F")), str(prows))

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("skip params hold in both directions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
