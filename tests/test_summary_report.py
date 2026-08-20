#!/usr/bin/env python3
"""Tests for annotate_snps' run summary TSV (05_report/<prefix>.summary.tsv).

    python3 tests/test_summary_report.py

Two layers, cheapest first:

  script     bin/summary_report.py run directly on synthetic fixtures. This is where the
             numbers are checked: dedup of the (variant x carrier-sample) multianno rows,
             the consequence and ClinVar breakdowns, the threshold flags, sex from X/Y
             depth ratios, and the tolerant path -- a multianno with none of the optional
             columns (the exact header test_skip_params' fake table_annovar writes) must
             produce a report and exit 0, because a report must never kill the run.
  pipeline   the real DAG with test_skip_params' fake tools (imported, not copied), one
             completing run with QC + ancestry on. This is where the wiring is checked:
             summary_report actually runs, the somalier channels reach it through
             toList(), and the TSV lands in 05_report/ under the documented name.

Needs python3 only for the script layer; the pipeline layer is skipped without nextflow,
same as the rest of the harness.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
SCRIPT = ROOT / "nf" / "annotate_snps" / "bin" / "summary_report.py"

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run_script(tmp, out_name, extra):
    cmd = [sys.executable, str(SCRIPT),
           "--cohort", "TEST", "--data-type", "wgs",
           "--run-date", "20260627", "--clinvar-date", "clinvar_20260627",
           "--out", out_name] + extra
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp))
    rows = {}
    out = tmp / out_name
    if out.is_file():
        for line in out.read_text().splitlines()[1:]:
            section, metric, sample, value = line.split("\t")
            rows[(section, metric, sample)] = value
    return r, rows


# --- fixtures: two samples, three variants, one shared -------------------------------------

STATS = """\
SN\t0\tnumber of samples:\t2
SN\t0\tnumber of records:\t100
SN\t0\tnumber of SNPs:\t80
SN\t0\tnumber of indels:\t20
SN\t0\tnumber of multiallelic sites:\t3
# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)
TSTV\t0\t60\t20\t3.00\t60\t20\t3.00
# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing
PSC\t0\tSAMP_F\t50\t10\t40\t30\t10\t5\t31.2\t4\t0\t0\t1
PSC\t0\tSAMP_M\t55\t20\t25\t28\t12\t6\t29.8\t2\t0\t0\t0
"""

MULTIANNO_HEADER = ("Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tExonicFunc.refGene\tCLNSIG"
                    "\tCADDv1.6\tAF\tSpliceAI_DS_max\tOtherinfo1\tOtherinfo2\n")
# v1 in both samples (must count once), v2/v3 in one each. v1 missense, pathogenic,
# CADD 25, gnomAD-novel, splice-high; v2 synonymous, benign, CADD 10; v3 intronic,
# no exonic func.
MULTIANNO = MULTIANNO_HEADER + (
    "chr1\t100\t100\tA\tG\texonic\tnonsynonymous SNV\tPathogenic\t25.1\t.\t0.91\tx\tSAMP_F\n"
    "chr1\t100\t100\tA\tG\texonic\tnonsynonymous SNV\tPathogenic\t25.1\t.\t0.91\tx\tSAMP_M\n"
    "chr1\t200\t200\tC\tT\texonic\tsynonymous SNV\tBenign\t10.0\t0.005\t0.01\tx\tSAMP_F\n"
    "chr2\t300\t300\tG\tA\tintronic\t.\t.\t3.2\t0.002\t.\tx\tSAMP_M\n"
)

SAMPLES_TSV = ("#family_id\tsample_id\tdepth_mean\tX_depth_mean\tY_depth_mean\n"
               "FAM\tSAMP_M\t30\t15\t12\n"
               "FAM\tSAMP_F\t30\t30\t0.3\n")
PLOIDY = ("#sample\tx_copies\ty_copies\tped_sex\tagreement\n"
          "SAMP_M\t1\t1\t1\tAGREES\n"
          "SAMP_F\t2\t0\t1\tDISAGREES\n")
ANCESTRY = ("#sample_id\tpredicted_ancestry\tgiven_ancestry\tEUR_prob\tAFR_prob\n"
            "SAMP_M\tEUR\t\t0.98\t0.01\n"
            "SAMP_F\tAFR\t\t0.02\t0.95\n")
PAIRS = ("#sample_a\tsample_b\trelatedness\texpected_relatedness\n"
         "SAMP_M\tSAMP_F\t0.49\t-1.0\n")


def script_full(tmp):
    print("script -- full inputs")
    for name, text in [("stats.txt", STATS), ("multi.tsv", MULTIANNO),
                       ("samples.tsv", SAMPLES_TSV), ("ploidy.tsv", PLOIDY),
                       ("ancestry.tsv", ANCESTRY), ("pairs.tsv", PAIRS)]:
        (tmp / name).write_text(text)
    r, rows = run_script(tmp, "full.tsv", [
        "--stats", "stats.txt", "--multianno", "multi.tsv",
        "--samples-tsv", "samples.tsv", "--ploidy", "ploidy.tsv",
        "--ancestry", "ancestry.tsv", "--pairs", "pairs.tsv"])
    check("exits 0", r.returncode == 0, r.stderr[-400:])

    expect = {
        ("run", "cohort", "."): "TEST",
        ("run", "n_samples", "."): "2",
        ("cohort", "total_records", "."): "100",
        ("cohort", "snps", "."): "80",
        ("cohort", "indels", "."): "20",
        ("cohort", "multiallelic_sites", "."): "3",
        ("cohort", "ts_tv", "."): "3.00",
        ("sample", "n_het", "SAMP_F"): "40",
        ("sample", "het_hom_ratio", "SAMP_F"): "4.0000",
        ("sample", "mean_depth", "SAMP_M"): "29.8",
        ("sample", "n_singletons", "SAMP_M"): "2",
        # 3 unique variants across 4 rows -- the shared one counted once.
        ("rare_subset", "rare_variants", "."): "3",
        ("rare_subset", "rare_variant_sample_rows", "."): "4",
        ("rare_subset", "func:exonic", "."): "2",
        ("rare_subset", "func:intronic", "."): "1",
        ("rare_subset", "exonic_func:nonsynonymous_SNV", "."): "1",
        ("rare_subset", "exonic_func:synonymous_SNV", "."): "1",
        # Per-sample carrier counts via the content-recognized Otherinfo2 column.
        ("rare_subset", "rare_variants", "SAMP_F"): "2",
        ("rare_subset", "rare_variants", "SAMP_M"): "2",
        ("flags", "clnsig:Pathogenic", "."): "1",
        ("flags", "clnsig:Benign", "."): "1",
        ("flags", "clinvar_path_or_likely_path", "."): "1",
        ("flags", "cadd_ge_20", "."): "1",
        ("flags", "gnomad_novel", "."): "1",
        ("flags", "spliceai_ge_0.5", "."): "1",
        ("qc", "inferred_sex", "SAMP_M"): "male",
        ("qc", "inferred_sex", "SAMP_F"): "female",
        ("qc", "ped_sex", "SAMP_F"): "1",
        ("qc", "sex_agreement", "SAMP_F"): "DISAGREES",
        ("qc", "predicted_ancestry", "SAMP_M"): "EUR",
        ("qc", "ancestry_prob", "SAMP_M"): "0.98",
        ("qc", "predicted_ancestry", "SAMP_F"): "AFR",
        ("qc", "max_relatedness", "."): "0.49",
        ("qc", "max_relatedness_pair", "."): "SAMP_M,SAMP_F",
    }
    for key, want in expect.items():
        check("/".join(key), rows.get(key) == want, f"got {rows.get(key)!r}")
    check("no squirls flag row -- no squirls column in the header",
          not any(m.startswith("squirls") for _s, m, _x in rows), str(sorted(rows)))


def script_degraded(tmp):
    print("script -- fake-harness multianno header, empty stats, no somalier files")
    # The exact header test_skip_params' fake table_annovar writes. Every optional column
    # and all somalier inputs are absent -- the report must still be written, and exit 0.
    (tmp / "bare.tsv").write_text("Chr\tStart\tEnd\tRef\tAlt\tCADD_phred\n"
                                  "chr1\t100\t100\tA\tG\t9.566\n")
    (tmp / "empty_stats.txt").write_text("")
    r, rows = run_script(tmp, "degraded.tsv", [
        "--stats", "empty_stats.txt", "--multianno", "bare.tsv"])
    check("exits 0", r.returncode == 0, r.stderr[-400:])
    check("run rows still present", rows.get(("run", "cohort", ".")) == "TEST", str(rows))
    check("the variant is still counted",
          rows.get(("rare_subset", "rare_variants", ".")) == "1", str(rows))
    check("CADD_phred is recognized as the CADD column",
          rows.get(("flags", "cadd_ge_20", ".")) == "0", str(rows))
    check("no cohort rows from empty stats",
          not any(s == "cohort" for s, _m, _x in rows), str(sorted(rows)))
    check("no clinvar flag rows without CLNSIG",
          not any(s == "flags" and (m.startswith("clnsig") or m.startswith("clinvar"))
                  for s, m, _x in rows), str(sorted(rows)))
    check("absences are noted on stderr, not fatal",
          "NOTE" in r.stderr, r.stderr[-300:])


def pipeline(tmp):
    print("pipeline -- summary_report runs and publishes to 05_report/")
    import test_skip_params as tsp
    binf = tsp.build_env(tmp)
    # QC + ancestry on, --skip_spliceai_squirls so the run COMPLETES locally (see
    # test_skip_params on why the unskipped DAG cannot finish outside the containers).
    _trace, _named, out = tsp.run_pipeline(
        tmp, binf,
        ["--skip_spliceai_squirls",
         "--somalier_sites", "sites.vcf.gz", "--ped", "cohort.ped",
         "--somalier_labels", "ancestry-labels.tsv",
         "--somalier_1kg_dir", "1kg"], "report")
    check("the run completes",
          "[FAILED]" not in out and "ERROR" not in out, f"...{out.strip()[-500:]}")

    tsv = tmp / "results_report" / "05_report" / "TEST_wgs_20260627.summary.tsv"
    check("05_report/TEST_wgs_20260627.summary.tsv is published", tsv.is_file(),
          str(sorted(p.name for p in tsv.parent.glob("*")) if tsv.parent.is_dir()
              else "no 05_report"))
    if not tsv.is_file():
        return
    rows = {}
    for line in tsv.read_text().splitlines()[1:]:
        section, metric, sample, value = line.split("\t")
        rows[(section, metric, sample)] = value
    check("header is the tidy four columns",
          tsv.read_text().startswith("#section\tmetric\tsample\tvalue\n"),
          tsv.read_text()[:80])
    # The fake somalier gives SAMP_M X ratio 0.5 / Y 0.4 and SAMP_F 1.0 / 0.01 -- the same
    # depths test_skip_params pins the ploidy table with.
    check("inferred sex reaches the report",
          rows.get(("qc", "inferred_sex", "SAMP_M")) == "male"
          and rows.get(("qc", "inferred_sex", "SAMP_F")) == "female", str(rows))
    check("the PED cross-check reaches the report",
          rows.get(("qc", "sex_agreement", "SAMP_F")) == "DISAGREES", str(rows))
    check("predicted ancestry reaches the report",
          rows.get(("qc", "predicted_ancestry", "SAMP_M")) == "EUR", str(rows))
    check("observed relatedness reaches the report",
          ("qc", "max_relatedness", ".") in rows, str(sorted(rows)))
    check("the multianno rare subset is counted",
          rows.get(("rare_subset", "rare_variants", ".")) == "1", str(rows))
    # The fake bcftools has no stats subcommand, so the cohort section must degrade away
    # rather than fail the task -- that IS the tolerant contract, exercised in-DAG.
    check("run metadata rows are present",
          rows.get(("run", "cohort", ".")) == "TEST", str(rows))


def main():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        script_full(tmp)
        script_degraded(tmp)
        if shutil.which("nextflow"):
            pipeline(tmp)
        else:
            print("SKIP: nextflow not on PATH -- pipeline layer not run")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("summary report holds at both layers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
