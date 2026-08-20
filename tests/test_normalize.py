#!/usr/bin/env python3
"""Regression tests for SV field normalization.

Covers nf/annotate_svs/assets/normalize_records.awk directly -- no Nextflow, no container.

These matter more than they look. Normalization used to be `svtk standardize`'s job, but
svtk cannot run (its bundled template is invalid VCF under modern pysam, in every current
biocontainer build), so this awk is the pipeline's ONLY standardization and it runs over
manta and delly too. Every failure here is a silent one: a wrong SVLEN sign or a missing
END does not error, it just stops the record matching anything, and the same deletion gets
reported twice from two callers.

    python tests/test_normalize.py

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

AWK = (Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
       / "normalize_records.awk")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Secondary chromosome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1"""

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def norm(info, alt="<DEL>", pos=1000):
    vcf = (HEADER + "\n"
           + f"chr1\t{pos}\trec\tN\t{alt}\t99\tPASS\t{info}\tGT\t0/1\n")
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.vcf"
        src.write_text(vcf)
        res = subprocess.run(["awk", "-f", str(AWK), str(src)],
                             capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"awk failed: {res.stderr.strip()[:300]}")
    rec = [l for l in res.stdout.splitlines() if not l.startswith("#")][0]
    return rec.split("\t")[7]


def val(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("SVLEN sign is normalized by type")
    # The failure mode: two callers disagree on the sign, size-similarity matching fails, and
    # one deletion becomes two clusters.
    check("a positive DEL SVLEN becomes negative",
          val(norm("SVTYPE=DEL;END=2000;SVLEN=1000"), "SVLEN") == "-1000")
    check("an already-negative DEL SVLEN is left negative",
          val(norm("SVTYPE=DEL;END=2000;SVLEN=-1000"), "SVLEN") == "-1000")
    check("a negative DUP SVLEN becomes positive",
          val(norm("SVTYPE=DUP;END=2000;SVLEN=-1000", alt="<DUP>"), "SVLEN") == "1000")
    check("INV SVLEN is positive",
          val(norm("SVTYPE=INV;END=2000;SVLEN=-1000", alt="<INV>"), "SVLEN") == "1000")

    print("END is resolved rather than left as a 1 bp event")
    check("a missing END is computed from SVLEN",
          val(norm("SVTYPE=DEL;SVLEN=-1000"), "END") == "2000",
          "without END every downstream tool treats this as 1 bp and it never merges")
    check("an END collapsed to POS is recomputed",
          val(norm("SVTYPE=DEL;END=1000;SVLEN=-1000"), "END") == "2000")
    check("a valid END is preserved",
          val(norm("SVTYPE=DEL;END=2500;SVLEN=-1500"), "END") == "2500")
    check("SVLEN is derived from END when only END is given",
          val(norm("SVTYPE=DEL;END=2000"), "SVLEN") == "-1000")

    print("SVTYPE is derived from a symbolic ALT when absent")
    check("<DEL> gives SVTYPE=DEL", val(norm("END=2000;SVLEN=-1000"), "SVTYPE") == "DEL")
    check("<DUP:TANDEM> gives SVTYPE=DUP",
          val(norm("END=2000;SVLEN=1000", alt="<DUP:TANDEM>"), "SVTYPE") == "DUP",
          "the subtype suffix must not leak into SVTYPE")
    check("an explicit SVTYPE is not overwritten by the ALT",
          val(norm("SVTYPE=INS;END=2000;SVLEN=1000", alt="<DUP>"), "SVTYPE") == "INS")

    print("INS records keep their END -- an insertion consumes no reference")
    # END == POS is spec-correct for an insertion, NOT a collapsed 1 bp event. Deriving
    # END = POS + SVLEN gave a 500 bp INS a 500 bp reference span it never occupied, which
    # then clustered against real DELs and DUPs at 0.95 reciprocal overlap and made AnnotSV
    # annotate genes the variant does not touch.
    ins = norm("SVTYPE=INS;END=1000;SVLEN=500", alt="<INS>", pos=1000)
    check("INS END is not pushed out by SVLEN", val(ins, "END") == "1000", f"INFO={ins}")
    check("INS SVLEN survives", val(ins, "SVLEN") == "500", f"INFO={ins}")
    check("INS SVLEN sign is still normalized positive",
          val(norm("SVTYPE=INS;END=1000;SVLEN=-500", alt="<INS>"), "SVLEN") == "500")

    print("symbolic ALTs containing digits are typed, not skipped")
    # <CN0>/<CN2>/<CN3> are what gCNV and GATK-SV emit for multi-allelic CNVs, and
    # <INS:ME:L1>/<INS:ME:SVA> are standard mobile-element calls. An alpha-only character
    # class declined to type every one of them, and an untyped record is exactly what the
    # type-aware mergers drop or mis-cluster.
    cn0 = norm("END=5000", alt="<CN0>", pos=3000)
    check("<CN0> is typed as CNV, not as the literal CN0",
          val(cn0, "SVTYPE") == "CNV",
          f"INFO={cn0} -- SVDB and AnnotSV both switch on SVTYPE, so a bogus value is "
          "worse than an absent one")
    check("and its SVLEN is then derived", val(cn0, "SVLEN") == "2000", f"INFO={cn0}")
    me = norm("END=1000;SVLEN=6000", alt="<INS:ME:L1>", pos=1000)
    check("<INS:ME:L1> gives SVTYPE=INS", val(me, "SVTYPE") == "INS", f"INFO={me}")
    check("and is exempt from END derivation like any other INS",
          val(me, "END") == "1000", f"INFO={me}")

    print("END before POS is reported rather than passed on silently")
    # Repairable when SVLEN is there; when it is not, htslib rejects the record later in a
    # message naming neither the caller nor the record, so say which record here.
    check("END < POS is repaired from SVLEN when possible",
          val(norm("SVTYPE=DEL;END=100;SVLEN=-1000", pos=5000), "END") == "6000")
    bad = subprocess.run(
        ["awk", "-f", str(AWK)],
        input="chr1\t5000\trec\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL;END=100\tGT\t0/1\n",
        capture_output=True, text=True)
    check("an unrepairable END < POS warns on stderr naming the record",
          "WARNING" in bad.stderr and "rec" in bad.stderr,
          f"stderr={bad.stderr.strip()[:200]!r}")
    check("and the record is still emitted, not dropped",
          len([l for l in bad.stdout.splitlines() if not l.startswith("#")]) == 1)

    print("BND records are left alone")
    # A breakend has no span. Computing END or SVLEN for one invents data, and dropping
    # them loses the translocation silently.
    bnd = norm("SVTYPE=BND;CHR2=chr2;END=500", alt="N[chr2:500[")
    check("BND END is untouched", val(bnd, "END") == "500", f"INFO={bnd}")
    check("BND gains no SVLEN", val(bnd, "SVLEN") is None, f"INFO={bnd}")
    check("BND CHR2 survives", val(bnd, "CHR2") == "chr2", f"INFO={bnd}")

    print("everything else is preserved")
    got = norm("SVTYPE=DEL;END=2000;SVLEN=1000;CHR2=chr1;ALGORITHMS=manta")
    check("unrelated INFO keys survive",
          val(got, "CHR2") == "chr1" and val(got, "ALGORITHMS") == "manta",
          f"INFO={got}")
    check("no key is duplicated",
          len([kv.split("=")[0] for kv in got.split(";")])
          == len({kv.split("=")[0] for kv in got.split(";")}), f"INFO={got}")

    got = norm(".", alt="<DEL>")
    check("an empty INFO still gains SVTYPE from the ALT",
          val(got, "SVTYPE") == "DEL", f"INFO={got!r}")

    got = norm("SVTYPE=DEL;END=2000;SVLEN=-1000")
    check("an already-correct record is unchanged",
          val(got, "SVTYPE") == "DEL" and val(got, "END") == "2000"
          and val(got, "SVLEN") == "-1000", f"INFO={got}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("normalization holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
