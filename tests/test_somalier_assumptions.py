#!/usr/bin/env python3
"""The somalier behaviours the ploidy design rests on, asserted rather than assumed.

    python3 tests/test_somalier_assumptions.py

Same purpose as tests/test_svdb_assumptions.py: `nf/shared/assets/check_somalier.awk` reads four
columns out of somalier's samples table BY NAME and derives every sample's chrX/chrY copy
number from them. If a somalier release renames a column or changes what it measures, the
depth filter on the sex chromosomes silently changes meaning -- so the contract is pinned
here against the real tool rather than against a hand-written fixture.

Runs somalier in Docker over a toy three-contig genome built here: chr1 at full depth,
chrX at half depth for the male and full for the female, chrY present only for the male.
That is the entire signal the ploidy call is made from, so a toy genome is enough to pin it.

Needs docker and samtools/bgzip/tabix (environment.yml, conda env `annotate-svs`).
"""

import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent.parent
AWK = ROOT / "nf" / "shared" / "assets" / "check_somalier.awk"
# The tag pinned in nf/shared/modules/qc/somalier.nf. Kept in step with it by hand;
# this test is what would notice the two diverging.
IMAGE = "quay.io/biocontainers/somalier:0.2.19--h0c29559_0"

CONTIGS = {"chr1": 40000, "chrX": 20000, "chrY": 10000}
READ_LEN = 100

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run(cmd, **kw):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(str(c) for c in cmd)}\n{r.stdout[-2000:]}\n"
                           f"{r.stderr[-2000:]}")
    return r


def build_inputs(tmp):
    random.seed(7)
    seqs = {c: "".join(random.choice("ACGT") for _ in range(n))
            for c, n in CONTIGS.items()}
    fa = tmp / "ref.fa"
    with open(fa, "w") as fh:
        for c, s in seqs.items():
            fh.write(f">{c}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    run(["samtools", "faidx", fa])

    # Sites every 500 bp. somalier genotypes these and reports mean depth over them per
    # contig group -- which is where X_depth_mean and Y_depth_mean come from.
    sites = []
    for c, n in CONTIGS.items():
        for pos in range(1000, n - 1000, 500):
            ref = seqs[c][pos - 1]
            sites.append((c, pos, ref, {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]))
    vcf = tmp / "sites.vcf"
    with open(vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, n in CONTIGS.items():
            fh.write(f"##contig=<ID={c},length={n}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, pos, ref, alt in sites:
            fh.write(f"{c}\t{pos}\t{c}_{pos}\t{ref}\t{alt}\t100\tPASS\t.\n")
    run(["bgzip", "-f", str(vcf)])
    run(["tabix", "-f", "-p", "vcf", str(vcf) + ".gz"])
    return seqs, sites


def write_bam(tmp, seqs, sites, sample, depths, alt_frac_x):
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines += [f"@SQ\tSN:{c}\tLN:{n}" for c, n in CONTIGS.items()]
    lines.append(f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA")
    body, n_read = [], 0
    for c, n in CONTIGS.items():
        if depths[c] == 0:
            continue
        step = max(1, READ_LEN // depths[c])
        for pos in range(1, n - READ_LEN, step):
            seq = list(seqs[c][pos - 1:pos - 1 + READ_LEN])
            if c == "chrX" and alt_frac_x and random.random() < alt_frac_x:
                for sc, sp, _ref, alt in sites:
                    if sc == c and pos <= sp < pos + READ_LEN:
                        seq[sp - pos] = alt
            n_read += 1
            body.append((c, pos, "\t".join([
                f"r{n_read}", "0", c, str(pos), "60", f"{READ_LEN}M", "*", "0", "0",
                "".join(seq), "I" * READ_LEN, f"RG:Z:{sample}"])))
    order = {c: i for i, c in enumerate(CONTIGS)}
    body.sort(key=lambda t: (order[t[0]], t[1]))
    sam = tmp / f"{sample}.sam"
    sam.write_text("\n".join(lines + [b[2] for b in body]) + "\n")
    run(["samtools", "view", "-b", "-o", tmp / f"{sample}.bam", sam])
    run(["samtools", "index", tmp / f"{sample}.bam"])


def read_tsv(path):
    rows = path.read_text().splitlines()
    hdr = [h.lstrip("#") for h in rows[0].split("\t")]
    return hdr, [dict(zip(hdr, r.split("\t"))) for r in rows[1:] if r.strip()]


def ancestry_output_naming():
    """somalier ancestry's output prefix, which decides the filenames the pipeline renames.

    `somalier ancestry` appends a FIXED `.somalier-ancestry.<ext>` to its `-o` prefix, and that
    prefix already defaults to `somalier-ancestry` -- so with no `-o` it writes
    `somalier-ancestry.somalier-ancestry.tsv`. somalier_ancestry renames those files, and
    tests/test_skip_params.py's fake reproduces the doubled names, so both go stale silently if
    the default ever changes. That defect shipped once already: the process declared output
    files somalier never wrote, and every ancestry run died on "Missing output file(s)".
    """
    r = subprocess.run(["docker", "run", "--rm", IMAGE, "somalier", "ancestry", "--help"],
                       capture_output=True, text=True)
    text = r.stdout + r.stderr
    check("ancestry still takes -o/--output-prefix",
          "--output-prefix" in text, text[-300:])
    check("and it still defaults to somalier-ancestry",
          "default: somalier-ancestry" in text,
          "if this moved, update the mv in nf/shared/modules/qc/somalier.nf AND the fake in "
          f"tests/test_skip_params.py: {text[-300:]}")


def main():
    for tool in ("docker", "samtools", "bgzip", "tabix"):
        if not shutil.which(tool):
            print(f"SKIP: {tool} not available")
            return 0

    ancestry_output_naming()

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        seqs, sites = build_inputs(tmp)
        # The male's chrX runs at half the autosomal depth and his chrY is present; the
        # female's chrX matches chr1 and she has no chrY reads at all. Everything below is
        # derived from exactly that.
        write_bam(tmp, seqs, sites, "MALE1", {"chr1": 20, "chrX": 10, "chrY": 8}, 0.0)
        write_bam(tmp, seqs, sites, "FEM1", {"chr1": 20, "chrX": 20, "chrY": 0}, 0.5)
        (tmp / "cohort.ped").write_text(
            "F1\tMALE1\t0\t0\t1\t2\nF1\tFEM1\t0\t0\t2\t1\n")

        docker = ["docker", "run", "--rm", "-v", f"{tmp}:/d", "-w", "/d", IMAGE]
        for s in ("MALE1", "FEM1"):
            run(docker + ["somalier", "extract", "-d", "/d", "--sites", "sites.vcf.gz",
                          "-f", "ref.fa", f"{s}.bam"])

        print("somalier extract")
        check("names its output after the read-group sample, not the file",
              (tmp / "MALE1.somalier").is_file() and (tmp / "FEM1.somalier").is_file(),
              "the module asserts this filename; if it changes, that guard misfires")

        run(docker + ["bash", "-lc", "somalier relate --ped cohort.ped -o cohort "
                                     "MALE1.somalier FEM1.somalier"])

        print("samples.tsv carries the columns check_somalier.awk reads by name")
        hdr, rows = read_tsv(tmp / "cohort.samples.tsv")
        for col in ("sample_id", "depth_mean", "X_depth_mean", "Y_depth_mean"):
            check(f"{col} is present", col in hdr, str(hdr))
        by = {r["sample_id"]: r for r in rows}

        print("the depth ratios separate the two karyotypes")
        m, f = by["MALE1"], by["FEM1"]
        mx = float(m["X_depth_mean"]) / float(m["depth_mean"])
        fx = float(f["X_depth_mean"]) / float(f["depth_mean"])
        my = float(m["Y_depth_mean"]) / float(m["depth_mean"])
        fy = float(f["Y_depth_mean"]) / float(f["depth_mean"])
        check("one X reads at about half the autosomal depth", mx <= 0.65, f"{mx:.2f}")
        check("two X read at about the autosomal depth", fx >= 0.80, f"{fx:.2f}")
        check("a Y that is present clears the cut point", my >= 0.15, f"{my:.2f}")
        check("a Y that is absent falls under it", fy <= 0.05, f"{fy:.2f}")

        print("pairs.tsv gives an expectation to check relatedness against")
        phdr, prows = read_tsv(tmp / "cohort.pairs.tsv")
        for col in ("sample_a", "sample_b", "relatedness", "expected_relatedness"):
            check(f"{col} is present", col in phdr, str(phdr))
        check("an unrelated-by-PED pair is marked -1, not 0",
              bool(prows) and float(prows[0]["expected_relatedness"]) < 0,
              str(prows[:1]) +
              " -- the awk skips negatives; a 0 would make every unrelated pair a warning")

        print("check_somalier.awk over the real output")
        run(["awk",
             "-v", f"SAMPLES={tmp}/cohort.samples.tsv",
             "-v", f"PAIRS={tmp}/cohort.pairs.tsv",
             "-v", f"PED={tmp}/cohort.ped",
             "-v", f"OUT_PLOIDY={tmp}/ploidy.tsv",
             "-v", f"OUT_PED={tmp}/cohort.inferred.ped",
             "-f", str(AWK)], stdin=subprocess.DEVNULL)
        ploidy = {r[0]: r[1:] for r in
                  (l.split("\t") for l in (tmp / "ploidy.tsv").read_text().splitlines()
                   if not l.startswith("#") and l.strip())}
        check("the male comes out 1 X, 1 Y", ploidy["MALE1"][:2] == ["1", "1"],
              str(ploidy))
        check("the female comes out 2 X, 0 Y", ploidy["FEM1"][:2] == ["2", "0"],
              str(ploidy))
        check("both agree with the PED", ploidy["MALE1"][3] == "AGREES"
              and ploidy["FEM1"][3] == "AGREES", str(ploidy))
        check("the inferred PED keeps the sexes it measured",
              [l.split("\t")[4] for l in
               (tmp / "cohort.inferred.ped").read_text().splitlines()] == ["1", "2"],
              (tmp / "cohort.inferred.ped").read_text())

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("somalier assumptions hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
