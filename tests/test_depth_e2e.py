#!/usr/bin/env python3
"""End-to-end regression test for the sample-to-alignment pairing in Phase 3.

    python tests/test_depth_e2e.py

Unlike test_depth_merge.py, which exercises merge_depth.awk on its own, this one runs the
actual pipeline with --alignments. The defect it guards against lives in Nextflow channel
wiring and cannot be reproduced any other way.

WHY THIS EXISTS. depth_svs splits the cohort VCF per sample, hands each single-sample VCF to
duphold along with that sample's alignment, and joins the results back. Until this test there
was no automated coverage of that path at all, and a rewiring of it shipped broken: the sites
VCF and the alignment reached duphold down two separate routes with nothing checking they
described the same sample. Depth attached to the wrong sample does not error. It produces a
complete, plausible, wrong answer, and then drives the depth filter from it -- DEL_DHFFC and
DUP_DHBFC are depth thresholds, so a swap turns a real deletion into a filtered artefact and
an artefact into a confident call.

HOW IT IS DETECTED. The two samples' alignments are deliberately not interchangeable. S1 has
a real deletion over chr1:1000-2000 -- no reads there, normal flanks -- and S2 has uniform
coverage throughout. duphold's DHFFC is depth over the SV relative to its flanks, so the one
record in the cohort VCF must come out with a LOW DHFFC for S1 and a normal one for S2. Swap
the two BAMs and both assertions invert. A test that merely checked "every sample got a
number" would pass on the swap, which is the failure this is for.

ALSO COVERED. A sheet whose `sample` column matches nothing in the cohort. That used to
continue: the depth stages ran over nothing and qc_somalier handed somalier an empty file
list, which 0.2.19 dies on with an IndexDefect naming neither the sheet nor the cohort.

Needs samtools (environment.yml, conda env `annotate-svs`) to build the BAMs.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
MAIN = ROOT / "nf" / "annotate_svs" / "main.nf"

CONTIG_LEN = 5000
READ_LEN = 100
# The deletion S1 carries and S2 does not. Reads are omitted here for S1 only.
DEL_START, DEL_END = 1000, 2000

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)


def write_reference(tmp):
    seq = ("ACGT" * (CONTIG_LEN // 4))[:CONTIG_LEN]
    fa = tmp / "ref.fa"
    fa.write_text(">chr1\n"
                  + "\n".join(seq[i:i + 60] for i in range(0, CONTIG_LEN, 60)) + "\n")
    run(["samtools", "faidx", fa])
    return fa, seq


def write_bam(tmp, seq, sample, deleted):
    """Coordinate-sorted by construction: positions ascend, so no sort step is needed."""
    lines = ["@HD\tVN:1.6\tSO:coordinate",
             f"@SQ\tSN:chr1\tLN:{CONTIG_LEN}",
             f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"]
    n = 0
    for pos in range(1, CONTIG_LEN - READ_LEN, 20):
        if deleted and DEL_START <= pos <= DEL_END:
            continue
        n += 1
        lines.append("\t".join([
            f"r{n}", "0", "chr1", str(pos), "60", f"{READ_LEN}M", "*", "0", "0",
            seq[pos - 1:pos - 1 + READ_LEN], "I" * READ_LEN, f"RG:Z:{sample}"]))
    sam = tmp / f"{sample}.sam"
    sam.write_text("\n".join(lines) + "\n")
    bam = tmp / f"{sample}.bam"
    run(["samtools", "view", "-b", "-o", bam, sam])
    run(["samtools", "index", bam])
    return bam


def write_vcf(tmp, samples, caller="manta"):
    """One joint DEL over the window S1 is missing reads in.

    Written once per caller. Two callers over the same sample set is what makes axis A merge
    them, and an SVDB merge is what puts the per-input blobs into the cohort INFO -- which is
    what split_sample strips and what this test has to actually exercise. A single-input
    fixture never merges, so it would only ever take the "no blobs" branch of that guard.
    """
    header = ("##fileformat=VCFv4.2\n"
              f"##contig=<ID=chr1,length={CONTIG_LEN}>\n"
              "##ALT=<ID=DEL,Description=\"Deletion\">\n"
              "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n"
              "##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n"
              "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"l\">\n"
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n")
    plain = tmp / f"{caller}.vcf"
    plain.write_text(
        header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
        + f"chr1\t{DEL_START}\t.\tN\t<DEL>\t99\tPASS\t"
          f"SVTYPE=DEL;END={DEL_END};SVLEN=-{DEL_END - DEL_START}\tGT\t0/1\t0/0\n")
    gz = tmp / f"{caller}.vcf.gz"
    gz.write_bytes(subprocess.run(["bgzip", "-c", str(plain)], capture_output=True).stdout)
    run(["tabix", "-p", "vcf", gz])
    return gz


def main():
    for tool in ("nextflow", "samtools", "bcftools", "bgzip", "tabix", "docker"):
        if not shutil.which(tool):
            print(f"SKIP: {tool} not available")
            return 0

    env = dict(os.environ, NXF_SYNTAX_PARSER="v2", NXF_ANSI_LOG="false")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        fa, seq = write_reference(tmp)
        # S1 carries the deletion, S2 does not. This asymmetry is the whole test.
        write_bam(tmp, seq, "S1", deleted=True)
        write_bam(tmp, seq, "S2", deleted=False)
        # Two callers over the same samples, so axis A merges and the cohort INFO carries
        # SVDB's per-input blobs -- the ones split_sample strips.
        rows = "".join(
            f"COH,{caller},true,{gz},{gz}.tbi\n"
            for caller, gz in ((c, write_vcf(tmp, ["S1", "S2"], c))
                               for c in ("manta", "smoove")))
        (tmp / "vcfs.csv").write_text("sample_set,caller,joint,vcf,tbi\n" + rows)
        (tmp / "alignments.csv").write_text(
            "sample,alignment,alignment_index\n"
            f"S1,{tmp}/S1.bam,{tmp}/S1.bam.bai\n"
            f"S2,{tmp}/S2.bam,{tmp}/S2.bam.bai\n")
        (tmp / "cohort.ped").write_text("FAM\tS1\t0\t0\t1\t0\nFAM\tS2\t0\t0\t2\t0\n")

        print("depth evidence lands on the sample it was computed from")
        r = run(["nextflow", "-q", "run", MAIN, "-with-docker",
                 "--vcfs", tmp / "vcfs.csv", "--ped", tmp / "cohort.ped",
                 "--alignments", tmp / "alignments.csv",
                 "--alignment_reference", fa,
                 "--alignment_reference_index", f"{fa}.fai",
                 "--outdir", tmp / "results"],
                cwd=str(tmp), env=env)
        errors = [l for l in (r.stdout + r.stderr).splitlines() if l.startswith("[ERROR]")]
        check("the pipeline runs clean", not errors, "; ".join(errors[:2]))
        if errors:
            return 1

        cohort = tmp / "results" / "05_filter" / "cohort.tagged.vcf.gz"
        if not cohort.is_file():
            check("a cohort VCF was produced", False, str(cohort))
            return 1

        q = run(["bcftools", "query", "-f", "[%SAMPLE=%DHFFC\t]\n", cohort])
        fields = dict(f.split("=", 1) for f in q.stdout.strip().split("\t") if "=" in f)
        check("both samples carry a depth value", set(fields) == {"S1", "S2"}, str(fields))
        if set(fields) != {"S1", "S2"}:
            return 1
        check("neither sample's depth is missing",
              "." not in (fields["S1"], fields["S2"]), str(fields))
        if "." in (fields["S1"], fields["S2"]):
            return 1

        s1, s2 = float(fields["S1"]), float(fields["S2"])
        # The pairing assertion. Both halves are needed: `s1 < s2` alone would also hold if
        # every value collapsed toward zero, and the absolute bounds alone would not catch a
        # swap that kept both inside them.
        check("the sample WITH the deletion reads as depleted", s1 < 0.5, f"DHFFC S1={s1}")
        check("the sample WITHOUT it reads as normal", s2 > 0.5, f"DHFFC S2={s2}")
        check("and the two are not interchangeable", s1 < s2, f"S1={s1} S2={s2}")

        # split_sample strips SVDB's per-input blobs before duphold, because carried through
        # they make each uncompressed <sample>.dh.vcf tens of GB and the work dir multi-TB.
        # That strip must not reach the published VCF: the per-input provenance -- which
        # input supplied which genotype -- is answerable from this file and nowhere else.
        # Depth is joined back onto the full cohort VCF, so the blobs should all still be here.
        hdr = run(["bcftools", "view", "-h", cohort]).stdout
        blobs = [l.split("ID=")[1].split(",")[0] for l in hdr.splitlines()
                 if l.startswith("##INFO=<ID=")
                 and l.split("ID=")[1].split(",")[0].endswith(
                     ("_INFO", "_SAMPLE", "_FILTERS", "_CHROM", "_POS", "_QUAL", "_FORMAT"))]
        # If this fails the fixture stopped merging, and the strip branch is going untested.
        check("the fixture actually produced blobs to strip", blobs, "no <tag>_* INFO keys")
        check("and the published cohort VCF still carries them",
              bool(blobs) and "DHFFC" in hdr, f"blobs={blobs[:3]}")

        # Same fixtures, same BAMs -- only the sheet's sample column moves. A total mismatch
        # used to be treated like the partial one (reported at preflight, rows ignored), which
        # left in_cohort empty: depth ran over nothing and somalier was handed an empty file
        # list. Asserted on the message text, not just the exit code, because the whole point
        # of the change is that the operator can see WHICH side is wrong.
        print("an alignments sheet matching no cohort sample fails by name")
        (tmp / "wrong.csv").write_text(
            "sample,alignment,alignment_index\n"
            f"NOPE1,{tmp}/S1.bam,{tmp}/S1.bam.bai\n"
            f"NOPE2,{tmp}/S2.bam,{tmp}/S2.bam.bai\n")
        r = run(["nextflow", "-q", "run", MAIN, "-with-docker",
                 "--vcfs", tmp / "vcfs.csv", "--ped", tmp / "cohort.ped",
                 "--alignments", tmp / "wrong.csv",
                 "--alignment_reference", fa,
                 "--alignment_reference_index", f"{fa}.fai",
                 "--outdir", tmp / "results_wrong"],
                cwd=str(tmp), env=env)
        out = r.stdout + r.stderr
        check("the run fails", r.returncode != 0, f"rc={r.returncode}")
        check("it names both sides",
              "Sheet: NOPE1, NOPE2" in out and "Cohort: S1, S2" in out, out[-1500:])
        check("and no cohort VCF is published",
              not (tmp / "results_wrong" / "05_filter" / "cohort.tagged.vcf.gz").is_file())

    print()
    if failures:
        print(f"{len(failures)} FAILED: " + ", ".join(failures))
        return 1
    print("each sample's depth came from its own alignment")
    return 0


if __name__ == "__main__":
    sys.exit(main())
