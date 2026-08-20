#!/usr/bin/env python3
"""End-to-end test for the somalier QC stage and the ploidy table it feeds.

    python3 tests/test_ploidy_e2e.py

WHY THIS EXISTS. tests/test_somalier_assumptions.py pins what somalier reports, and
tests/test_filter_tags.py pins what the thresholds do with it. Neither runs the pipeline, so
neither would catch the wiring between them: qc_somalier reaching tag_filters at both of its
call sites, the ploidy table surviving into 07_qc, and the inferred PED reaching the Talos
tail. A defect there does not error -- the ploidy table would simply be absent or empty, and
the sex chromosomes would go back to being exempt without anything saying so.

WHAT IT ASSERTS. The wiring and the sex calls, not duphold's numbers: DHBFC is normalized
against GC-matched bins genome-wide, which on a 70 kb toy genome is not a quantity worth
asserting a threshold against. The threshold arithmetic is pinned by unit tests instead.

Needs nextflow, samtools, bcftools, bgzip, tabix and docker (conda env `annotate-svs`).
"""

import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
MAIN = ROOT / "nf" / "annotate_svs" / "main.nf"

CONTIGS = {"chr1": 40000, "chrX": 20000, "chrY": 10000}
READ_LEN = 100
# Outside GRCh38's PAR1 (which ends at 2,781,479) only in the sense that the record does not
# take the PAR branch; the toy contig is far shorter than the real chrX either way.
DUP_START, DUP_END = 4000, 6000

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)


def build_reference(tmp):
    random.seed(11)
    seqs = {c: "".join(random.choice("ACGT") for _ in range(n))
            for c, n in CONTIGS.items()}
    fa = tmp / "ref.fa"
    with open(fa, "w") as fh:
        for c, s in seqs.items():
            fh.write(f">{c}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    run(["samtools", "faidx", fa])
    return fa, seqs


def build_sites(tmp, seqs):
    vcf = tmp / "sites.vcf"
    with open(vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, n in CONTIGS.items():
            fh.write(f"##contig=<ID={c},length={n}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, n in CONTIGS.items():
            for pos in range(1000, n - 1000, 500):
                ref = seqs[c][pos - 1]
                alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
                fh.write(f"{c}\t{pos}\t{c}_{pos}\t{ref}\t{alt}\t100\tPASS\t.\n")
    run(["bgzip", "-f", str(vcf)])
    run(["tabix", "-f", "-p", "vcf", str(vcf) + ".gz"])
    return tmp / "sites.vcf.gz"


def write_bam(tmp, seqs, sample, depths):
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines += [f"@SQ\tSN:{c}\tLN:{n}" for c, n in CONTIGS.items()]
    lines.append(f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA")
    body, n_read = [], 0
    for c, n in CONTIGS.items():
        if depths[c] == 0:
            continue
        step = max(1, READ_LEN // depths[c])
        for pos in range(1, n - READ_LEN, step):
            n_read += 1
            body.append((c, pos, "\t".join([
                f"r{n_read}", "0", c, str(pos), "60", f"{READ_LEN}M", "*", "0", "0",
                seqs[c][pos - 1:pos - 1 + READ_LEN], "I" * READ_LEN, f"RG:Z:{sample}"])))
    order = {c: i for i, c in enumerate(CONTIGS)}
    body.sort(key=lambda t: (order[t[0]], t[1]))
    sam = tmp / f"{sample}.sam"
    sam.write_text("\n".join(lines + [b[2] for b in body]) + "\n")
    run(["samtools", "view", "-b", "-o", tmp / f"{sample}.bam", sam])
    run(["samtools", "index", tmp / f"{sample}.bam"])


def write_vcf(tmp, samples, caller):
    header = ["##fileformat=VCFv4.2"]
    header += [f"##contig=<ID={c},length={n}>" for c, n in CONTIGS.items()]
    header += ['##ALT=<ID=DUP,Description="Duplication">',
               '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">',
               '##INFO=<ID=END,Number=1,Type=Integer,Description="e">',
               '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="l">',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">']
    plain = tmp / f"{caller}.vcf"
    plain.write_text(
        "\n".join(header) + "\n"
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples) + "\n"
        + f"chrX\t{DUP_START}\t.\tN\t<DUP>\t99\tPASS\t"
          f"SVTYPE=DUP;END={DUP_END};SVLEN={DUP_END - DUP_START}\tGT\t0/1\t0/1\n")
    gz = tmp / f"{caller}.vcf.gz"
    gz.write_bytes(subprocess.run(["bgzip", "-c", str(plain)],
                                  capture_output=True).stdout)
    run(["tabix", "-p", "vcf", gz])
    return gz


def main():
    for tool in ("nextflow", "samtools", "bcftools", "bgzip", "tabix", "docker"):
        if not shutil.which(tool):
            print(f"SKIP: {tool} not available")
            return 0

    env = dict(os.environ, NXF_SYNTAX_PARSER="v2")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        fa, seqs = build_reference(tmp)
        sites = build_sites(tmp, seqs)
        # MALE1 has one X and a Y; FEM1 has two X and no Y. That asymmetry is the test.
        write_bam(tmp, seqs, "MALE1", {"chr1": 20, "chrX": 10, "chrY": 8})
        write_bam(tmp, seqs, "FEM1", {"chr1": 20, "chrX": 20, "chrY": 0})

        rows = ""
        for caller in ("manta", "smoove"):
            gz = write_vcf(tmp, ["MALE1", "FEM1"], caller)
            rows += f"COH,{caller},true,{gz},{gz}.tbi\n"
        (tmp / "vcfs.csv").write_text("sample_set,caller,joint,vcf,tbi\n" + rows)
        (tmp / "alignments.csv").write_text(
            "sample,alignment,alignment_index\n"
            f"MALE1,{tmp}/MALE1.bam,{tmp}/MALE1.bam.bai\n"
            f"FEM1,{tmp}/FEM1.bam,{tmp}/FEM1.bam.bai\n")
        # The PED calls FEM1 male. It is wrong on purpose: the data must win, and the
        # disagreement must be visible in the output rather than only in a log.
        (tmp / "cohort.ped").write_text(
            "FAM\tMALE1\t0\t0\t1\t0\nFAM\tFEM1\t0\t0\t1\t0\n")

        print("the somalier stage runs and its table reaches the filter")
        r = run(["nextflow", "-q", "run", MAIN, "-with-docker",
                 "--vcfs", tmp / "vcfs.csv", "--ped", tmp / "cohort.ped",
                 "--alignments", tmp / "alignments.csv",
                 "--alignment_reference", fa,
                 "--alignment_reference_index", f"{fa}.fai",
                 "--somalier_sites", sites,
                 "--outdir", tmp / "results"],
                cwd=str(tmp), env=env)
        errors = [l for l in (r.stdout + r.stderr).splitlines()
                  if l.startswith("[ERROR]")]
        check("the pipeline runs clean", not errors, "; ".join(errors[:2]))
        if errors:
            return 1

        qc = tmp / "results" / "07_qc"
        ploidy_file = qc / "ploidy.tsv"
        check("ploidy.tsv is published", ploidy_file.is_file(),
              str(sorted(p.name for p in qc.glob("*")) if qc.is_dir() else "no 07_qc"))
        if not ploidy_file.is_file():
            return 1

        ploidy = {row[0]: row[1:] for row in
                  (l.split("\t") for l in ploidy_file.read_text().splitlines()
                   if not l.startswith("#") and l.strip())}
        check("both samples are in it", set(ploidy) == {"MALE1", "FEM1"}, str(ploidy))
        check("the one-X sample is called hemizygous",
              ploidy.get("MALE1", [])[:2] == ["1", "1"], str(ploidy))
        check("the two-X sample is called diploid",
              ploidy.get("FEM1", [])[:2] == ["2", "0"], str(ploidy))
        check("the PED's wrong sex is reported as a disagreement",
              ploidy.get("FEM1", [])[3] == "DISAGREES", str(ploidy))
        check("and the sample the PED got right agrees",
              ploidy.get("MALE1", [])[3] == "AGREES", str(ploidy))

        print("the inferred PED carries the measured sex")
        inferred = qc / "cohort.inferred.ped"
        check("cohort.inferred.ped is published", inferred.is_file(), str(inferred))
        if inferred.is_file():
            lines = [l for l in inferred.read_text().splitlines() if l.strip()]
            sexes = {l.split("\t")[1]: l.split("\t")[4] for l in lines}
            check("the data overrode the PED", sexes.get("FEM1") == "2", str(sexes))
            check("family structure is untouched",
                  all(l.split("\t")[0] == "FAM" for l in lines), str(lines))

        print("relatedness QC is published beside it")
        for name in ("cohort.samples.tsv", "cohort.pairs.tsv", "cohort.html"):
            check(f"{name} is published", (qc / name).is_file(),
                  str(sorted(p.name for p in qc.glob("*"))))

        print("the ploidy used travels with the VCF")
        cohort = tmp / "results" / "05_filter" / "cohort.tagged.vcf.gz"
        check("a tagged cohort VCF was produced", cohort.is_file(), str(cohort))
        if cohort.is_file():
            hdr = run(["bcftools", "view", "-h", cohort]).stdout
            sample_sex = [l for l in hdr.splitlines() if l.startswith("##SAMPLE_SEX=")]
            check("one ##SAMPLE_SEX line per sample", len(sample_sex) == 2,
                  str(sample_sex))
            check("the disagreement is recorded in the VCF itself",
                  any("ID=FEM1" in l and "agreement=DISAGREES" in l for l in sample_sex),
                  str(sample_sex))

    print()
    if failures:
        print(f"{len(failures)} FAILED: " + ", ".join(failures))
        return 1
    print("ploidy reaches the filter, and the PED is checked rather than trusted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
