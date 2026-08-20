#!/usr/bin/env python3
"""Regression tests for the soft FILTER tags.

Covers nf/annotate_svs/assets/tag_filters.awk directly -- no Nextflow, no container.

The contract is "annotate hard, filter softly": every record is
kept, failing a criterion only writes a tag. A test suite is worth more here than almost
anywhere else in the pipeline, because a wrong threshold does not error -- it silently
changes which variants a reviewer ever sees.

The sex-chromosome case is the one to read first. A hemizygous male chrX sits near DHFFC
0.5, under the 0.7 deletion threshold, so a naive implementation calls every male chrX
event depth-supported. Until a PED supplies sample sex, depth tagging is autosome-only.

    python tests/test_filter_tags.py

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

AWK = (Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
       / "tag_filters.awk")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chrX,length=156040895>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=AF,Number=1,Type=Float,Description="Internal AF">
##INFO=<ID=NCALLER,Number=1,Type=Integer,Description="Callers">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DHFFC,Number=1,Type=Float,Description="flank fc">
##FORMAT=<ID=DHBFC,Number=1,Type=Float,Description="gc fc">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1"""

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def tag(info, fmt="GT", sample="0/1", chrom="chr1", pos=1000, ploidy=None, **kw):
    """Run the awk over a one-sample record; return (FILTER, full stdout).

    `ploidy` is a list of (sample, x_copies, y_copies, ped_sex, agreement) rows written to a
    temporary ploidy.tsv, or None for a run with no ploidy table at all -- the supported
    degraded configuration, not an error.
    """
    vcf = (HEADER + "\n"
           + f"{chrom}\t{pos}\trec_1\tN\t<DEL>\t99\t.\t{info}\t{fmt}\t{sample}\n")
    args = ["awk"]
    with tempfile.TemporaryDirectory() as td:
        if ploidy is not None:
            p = Path(td) / "ploidy.tsv"
            p.write_text("#sample\tx_copies\ty_copies\tped_sex\tagreement\n"
                         + "".join("\t".join(str(c) for c in row) + "\n" for row in ploidy))
            args += ["-v", f"PLOIDY={p}"]
        for k, v in kw.items():
            args += ["-v", f"{k}={v}"]
        args += ["-f", str(AWK)]
        src = Path(td) / "in.vcf"
        src.write_text(vcf)
        r = subprocess.run(args + [str(src)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"awk failed: {r.stderr.strip()[:300]}")
    rec = [l for l in r.stdout.splitlines() if not l.startswith("#")][0]
    return rec.split("\t")[6], r.stdout


def tag_cohort(info, genotypes, chrom="chr1", fmt="GT", **kw):
    """Run the awk over a cohort of len(genotypes) samples.

    Internal AF is computed from the GT columns, not read from INFO/AF, so exercising
    COMMON_INTERNAL needs a real cohort. A single-sample VCF cannot express any frequency
    below 0.5, which is why the awk declines to write the tag at that size at all.
    """
    names = "\t".join(f"S{i}" for i in range(1, len(genotypes) + 1))
    header = HEADER.rsplit("\tSAMP1", 1)[0] + "\t" + names
    vcf = (header + "\n"
           + f"{chrom}\t1000\trec_1\tN\t<DEL>\t99\t.\t{info}\t{fmt}\t"
           + "\t".join(genotypes) + "\n")
    args = ["awk"]
    for k, v in kw.items():
        args += ["-v", f"{k}={v}"]
    args += ["-f", str(AWK)]
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.vcf"
        src.write_text(vcf)
        r = subprocess.run(args + [str(src)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"awk failed: {r.stderr.strip()[:300]}")
    rec = [l for l in r.stdout.splitlines() if not l.startswith("#")][0]
    return rec.split("\t")[6], rec.split("\t")[7]


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("nothing failing means PASS, and every record survives")
    f, out = tag("SVTYPE=DEL;AF=0.001;NCALLER=3")
    check("a clean record is PASS", f == "PASS", f"got {f}")
    check("the record is still there",
          len([l for l in out.splitlines() if not l.startswith("#")]) == 1)

    print("population frequency")
    check("gnomAD SV AF at the threshold is tagged",
          tag("SVTYPE=DEL;gnomad_sv_AF=0.01;NCALLER=3")[0] == "COMMON_GNOMAD")
    check("gnomAD SV AF below it is not",
          tag("SVTYPE=DEL;gnomad_sv_AF=0.009;NCALLER=3")[0] == "PASS")
    check("the CNV callset is checked too, not just SV",
          tag("SVTYPE=DEL;gnomad_cnv_AF=0.5;NCALLER=3")[0] == "COMMON_GNOMAD",
          "two separate gnomAD files; only checking one silently misses half")
    check("a missing gnomAD AF is not treated as zero or as common",
          tag("SVTYPE=DEL;gnomad_sv_AF=.;NCALLER=3")[0] == "PASS")

    print("internal cohort frequency -- counted over ALL samples, not over called alleles")
    # 25 samples: the floor is 1/50 = 0.02, just under the 0.03 default, so the tag is usable
    # and a singleton is NOT common. Nine hets out of 25 is 9/50 = 0.18, which is.
    private = ["0/1"] + ["./."] * 24
    common = ["0/1"] * 9 + ["0/0"] * 16

    f, info = tag_cohort("SVTYPE=DEL;NCALLER=3", private)
    check("a private variant is not COMMON_INTERNAL, whatever fill-tags said",
          f == "PASS", f"got {f}",)
    check("and its computed frequency is written back as INTERNAL_AF",
          "INTERNAL_AF=0.02" in info, f"got {info}")

    check("a genuinely common variant is tagged",
          tag_cohort("SVTYPE=DEL;NCALLER=3", common)[0] == "COMMON_INTERNAL")
    # A doubleton, not the singleton: a singleton sits exactly AT the cohort floor, so no
    # threshold can tag it without also tagging every other carrier -- which is the condition
    # the awk refuses to write the tag under at all. Two hets in 25 is 2/50 = 0.04, comfortably
    # above the 0.02 floor, so the threshold genuinely decides the outcome here.
    doubleton = ["0/1", "0/1"] + ["./."] * 23
    check("the threshold is configurable -- below it, tagged",
          tag_cohort("SVTYPE=DEL;NCALLER=3", doubleton, INT_AF=0.035)[0] == "COMMON_INTERNAL")
    check("the threshold is configurable -- above it, not",
          tag_cohort("SVTYPE=DEL;NCALLER=3", doubleton, INT_AF=0.05)[0] == "PASS")

    # The regression this whole rewrite exists for. bcftools +fill-tags computes AN over
    # CALLED alleles only, and SVDB leaves a non-calling sample as "./.", so INFO/AF for a
    # private variant is 0.5 regardless of cohort size. Reading it tagged almost every
    # variant in the callset as common.
    check("INFO/AF is ignored even when it says 0.5",
          tag_cohort("SVTYPE=DEL;AF=0.5;NCALLER=3", private)[0] == "PASS",
          "reading INFO/AF here tagged every private variant COMMON_INTERNAL")

    # Below the size where the tag can discriminate, it is withheld rather than saturated.
    check("a cohort too small to express the threshold gets no tag at all",
          tag_cohort("SVTYPE=DEL;NCALLER=3", ["0/1", "0/0"])[0] == "PASS",
          "floor for 2 samples is 0.25, so every carrier would be 'common'")
    check("but the measured frequency is still reported",
          "INTERNAL_AF=0.25" in tag_cohort("SVTYPE=DEL;NCALLER=3", ["0/1", "0/0"])[1])

    print("multi-allelic records are judged on their worst allele")
    check("a Number=A gnomAD AF is not truncated to the first allele",
          tag("SVTYPE=DEL;gnomad_sv_AF=0.001,0.9;NCALLER=3")[0] == "COMMON_GNOMAD",
          "bare +0 reads 0.001 and lets a variant common at 0.9 through untagged")

    print("caller corroboration -- flagged, never dropped")
    check("a single-caller record is tagged",
          tag("SVTYPE=DEL;NCALLER=1")[0] == "NO_CALLER_SUPPORT")
    check("two callers is enough by default",
          tag("SVTYPE=DEL;NCALLER=2")[0] == "PASS")
    check("the record is kept regardless",
          len([l for l in tag("SVTYPE=DEL;NCALLER=1")[1].splitlines()
               if not l.startswith("#")]) == 1,
          "absence of support from a caller with a different sensitivity profile is "
          "weak evidence of absence")

    print("depth")
    check("a DEL whose depth did not drop is tagged",
          tag("SVTYPE=DEL;NCALLER=3", "GT:DHFFC", "0/1:0.95")[0] == "DEPTH_UNSUPPORTED")
    check("a DEL whose depth dropped is not",
          tag("SVTYPE=DEL;NCALLER=3", "GT:DHFFC", "0/1:0.05")[0] == "PASS")
    check("a DUP whose depth rose is not tagged",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.9")[0] == "PASS")
    check("a DUP whose depth did not rise is tagged",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0")[0] == "DEPTH_UNSUPPORTED")
    check("no depth measured means no opinion",
          tag("SVTYPE=DEL;NCALLER=3", "GT:DHFFC", "0/1:.")[0] == "PASS",
          "a sample with no alignment must not read as a failed depth check")
    check("a record with no depth fields at all is untouched",
          tag("SVTYPE=DEL;NCALLER=3")[0] == "PASS")

    print("sex chromosomes with no ploidy table are exempt, as they always were")
    check("a chrX record is not depth-tagged",
          tag("SVTYPE=DEL;NCALLER=3", "GT:DHFFC", "0/1:0.95", chrom="chrX")[0] == "PASS",
          "tagging chrX on autosomal thresholds is confidently wrong, not merely limited")
    check("chrY likewise",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="chrY")[0] == "PASS")
    check("chrX still gets the frequency tags, which do not depend on ploidy",
          tag_cohort("SVTYPE=DEL;NCALLER=3", ["0/1"] * 9 + ["0/0"] * 16,
                     chrom="chrX")[0] == "COMMON_INTERNAL")

    print("ploidy-aware sex chromosomes")
    male = [("SAMP1", 1, 1, 1, "AGREES")]
    female = [("SAMP1", 2, 0, 2, "AGREES")]
    unknown = [("SAMP1", ".", ".", 1, "UNKNOWN")]
    B = {"BUILD": "GRCh38"}

    # The defect this whole change exists for. A male chrX duplication of one copy into two
    # reads DHBFC ~1.0, because DHBFC is normalized genome-wide against a hemizygous
    # baseline of ~0.5. Against the flat 1.3 threshold every true male chrX DUP was tagged.
    check("a real male chrX DUP is no longer tagged",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="chrX",
              ploidy=male, **B)[0] == "PASS",
          "1.0 against a 0.5 baseline is a doubling; the scaled threshold is 0.65")
    check("but a male chrX DUP with no depth rise still is",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:0.5", chrom="chrX",
              ploidy=male, **B)[0] == "DEPTH_UNSUPPORTED",
          "the filter has to still work there, or scaling it just switched it off")
    check("the same value on a female chrX is tagged",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="chrX",
              ploidy=female, **B)[0] == "DEPTH_UNSUPPORTED",
          "two X copies means the autosomal threshold, and 1.0 is not a duplication")

    # DHFFC is flank-normalized and the flanks are equally hemizygous, so this threshold is
    # NOT scaled -- a male chrX deletion still reads near 0 and a non-deletion near 1.0.
    check("a male chrX DEL that dropped is supported",
          tag("SVTYPE=DEL;NCALLER=3", "GT:DHFFC", "0/1:0.05", chrom="chrX",
              ploidy=male, **B)[0] == "PASS")
    check("a male chrX DEL that did not drop is tagged",
          tag("SVTYPE=DEL;NCALLER=3", "GT:DHFFC", "0/1:0.95", chrom="chrX",
              ploidy=male, **B)[0] == "DEPTH_UNSUPPORTED")

    check("chrY in a sample with no Y is skipped, not judged",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:0.5", chrom="chrY",
              ploidy=female, **B)[0] == "PASS",
          "zero expected copies gives no threshold that means anything")
    check("chrY in a sample with a Y is judged on the hemizygous threshold",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="chrY",
              ploidy=male, **B)[0] == "PASS")

    check("an undetermined karyotype keeps the exemption",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:0.5", chrom="chrX",
              ploidy=unknown, **B)[0] == "PASS")
    check("a sample missing from the table keeps it too",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:0.5", chrom="chrX",
              ploidy=[("OTHER", 1, 1, 1, "AGREES")], **B)[0] == "PASS")
    check("autosomes are unaffected by the table",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", ploidy=male, **B)[0]
          == "DEPTH_UNSUPPORTED")

    print("PAR is diploid for everyone")
    # chrX:2,000,000 is inside GRCh38 PAR1 (10001-2781479), where a male carries two copies.
    check("a male PAR DUP is judged on the unscaled threshold",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="chrX", pos=2000000,
              ploidy=male, **B)[0] == "DEPTH_UNSUPPORTED",
          "PAR is diploid, so 1.0 is a baseline rather than a doubling")
    check("just outside PAR1 the hemizygous threshold applies again",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="chrX", pos=2800000,
              ploidy=male, **B)[0] == "PASS")
    check("PAR2 is recognised as well",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="chrX", pos=155800000,
              ploidy=male, **B)[0] == "DEPTH_UNSUPPORTED")
    check("a bare (non-chr) contig name is matched too",
          tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="X", pos=2000000,
              ploidy=male, **B)[0] == "DEPTH_UNSUPPORTED")
    # GRCh38 is the only build this pipeline supports, so every other one -- including a
    # coordinate-compatible-looking GRCh37 -- takes the fallback rather than PAR coordinates
    # that would be silently wrong for it.
    for build in ("GRCh37", "T2T-CHM13"):
        check(f"{build} falls back to exempting the sex chromosomes",
              tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:0.5", chrom="chrX",
                  ploidy=male, BUILD=build)[0] == "PASS",
              "PAR coordinates from the wrong build would be silent; none must not be")

    print("the ploidy used is recorded in the VCF")
    _, out = tag("SVTYPE=DUP;NCALLER=3", "GT:DHBFC", "0/1:1.0", chrom="chrX",
                 ploidy=[("SAMP1", 1, 1, 2, "DISAGREES")], **B)
    hdr = [l for l in out.splitlines() if l.startswith("##SAMPLE_SEX=")]
    check("a ##SAMPLE_SEX line is written per sample", len(hdr) == 1, str(hdr))
    check("it carries the copies actually used",
          "x_copies=1" in hdr[0] and "y_copies=1" in hdr[0], str(hdr))
    check("and a PED disagreement travels with the file",
          "agreement=DISAGREES" in hdr[0] and "ped_sex=2" in hdr[0], str(hdr),)
    check("no ploidy table means no such header",
          not [l for l in tag("SVTYPE=DEL;NCALLER=3")[1].splitlines()
               if l.startswith("##SAMPLE_SEX=")])

    print("BND and INS have no interval to measure")
    check("a BND is not depth-tagged",
          tag("SVTYPE=BND;NCALLER=3", "GT:DHFFC", "0/1:0.95")[0] == "PASS")
    check("an INS is not depth-tagged",
          tag("SVTYPE=INS;NCALLER=3", "GT:DHBFC", "0/1:1.0")[0] == "PASS")

    print("multiple failures accumulate")
    f, _ = tag_cohort("SVTYPE=DEL;gnomad_sv_AF=0.5;NCALLER=1",
                      ["0/1:0.95"] * 9 + ["0/0:1.0"] * 16, fmt="GT:DHFFC")
    parts = f.split(";")
    check("all four tags are present",
          sorted(parts) == sorted(["COMMON_GNOMAD", "COMMON_INTERNAL",
                                   "NO_CALLER_SUPPORT", "DEPTH_UNSUPPORTED"]),
          f"got {f}")
    check("and PASS is not among them", "PASS" not in parts, f"got {f}")

    print("headers")
    _, out = tag("SVTYPE=DEL;NCALLER=3")
    for t in ("COMMON_GNOMAD", "COMMON_INTERNAL", "DEPTH_UNSUPPORTED",
              "NO_CALLER_SUPPORT"):
        check(f"{t} is declared in the header",
              sum(1 for l in out.splitlines() if f"##FILTER=<ID={t}," in l) == 1)
    check("the declared threshold is the one applied",
          any("0.03" in l for l in out.splitlines()
              if "##FILTER=<ID=COMMON_INTERNAL" in l),
          "a header that documents a different number than the code uses is worse "
          "than no header")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("filter tagging holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
