#!/usr/bin/env python3
"""Regression tests for the Talos INFO/FORMAT keep-list.

Covers nf/annotate_svs/assets/talos_keep_info.awk directly -- no Nextflow, no container.

The keep-list exists because caller INFO baggage (PRPOS/PREND probability arrays,
caller-private keys) otherwise rides all the way into the 06_talos handoff file. The awk
prints the ^-prefixed KEEP form of `bcftools annotate -x`, restricted to tags the header
declares -- bcftools errors on a tag the file lacks. Its one dangerous failure mode is an
empty list: a bare "^" strips every INFO field, so the awk must exit nonzero with no
stdout instead.

    python tests/test_talos_strip.py

Needs `awk`; the bcftools round-trip section runs only when bcftools is on PATH (the test
env pins 1.24) and is skipped with a note otherwise.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ASSETS = Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
KEEP = ASSETS / "talos_keep_info.awk"

# What talos_schema.awk output looks like: required fields declared, plus the caller
# baggage this asset exists to remove. CHR2/END2 deliberately absent.
HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Len">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=AC,Number=1,Type=Integer,Description="AC">
##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">
##INFO=<ID=AF,Number=1,Type=Float,Description="AF">
##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Algs">
##INFO=<ID=SOFT_FILTERS,Number=.,Type=String,Description="Soft">
##INFO=<ID=gnomad_sv_AF,Number=1,Type=Float,Description="gnomAD">
##INFO=<ID=gnomad_v4.1_sv_AF,Number=1,Type=Float,Description="pop AF">
##INFO=<ID=PREDICTED_LOF,Number=.,Type=String,Description="LoF">
##INFO=<ID=PREDICTED_INTERGENIC,Number=0,Type=Flag,Description="Intergenic">
##INFO=<ID=PRPOS,Number=.,Type=String,Description="lumpy breakpoint probability">
##INFO=<ID=PREND,Number=.,Type=String,Description="lumpy breakpoint probability">
##INFO=<ID=SVDB_JUNK,Number=1,Type=String,Description="caller-private">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype likelihoods">
##FORMAT=<ID=PR,Number=.,Type=Integer,Description="manta paired-read support">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="lumpy split-read count -- the declaration manta's ref,alt pair violates">
##FORMAT=<ID=DHFFC,Number=1,Type=Float,Description="duphold">
##FORMAT=<ID=DHBFC,Number=1,Type=Float,Description="duphold">
##FORMAT=<ID=DHFC,Number=1,Type=Float,Description="duphold">
##FORMAT=<ID=DHBZ,Number=1,Type=Float,Description="duphold">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"""

# SR carries a manta-style two-valued "98,0" under lumpy's Number=1 declaration -- the
# exact shape hail refuses to load. The strip must remove it, header line and value.
RECORD = ("chr1\t100\tsv1\tN\t<DEL>\t.\tPASS\t"
          "SVTYPE=DEL;SVLEN=-500;END=600;AC=1;AN=2;AF=0.5;PREDICTED_LOF=GENE1;"
          "PRPOS=0.1,0.9;PREND=0.2,0.8;SVDB_JUNK=x;gnomad_sv_AF=0.001\t"
          "GT:PL:PR:SR:DHFFC:DHBFC:DHFC:DHBZ\t"
          "0/1:741,0,999:271,11:98,0:0.61:0.58:0.63:-1.2")

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run_keep(vcf_text, pop="gnomad_v4.1"):
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.vcf"
        src.write_text(vcf_text)
        return subprocess.run(["awk", "-v", f"GNOMAD_POP={pop}", "-f", str(KEEP), str(src)],
                              capture_output=True, text=True)


print("keep-list contents")
r = run_keep(HEADER + "\n" + RECORD + "\n")
check("awk exits 0", r.returncode == 0, r.stderr[:200])
out = r.stdout.strip()
check("one line, ^-prefixed", out.startswith("^") and "\n" not in out, out[:120])
tags = out.lstrip("^").split(",")
for want in ["INFO/SVTYPE", "INFO/SVLEN", "INFO/END", "INFO/AC", "INFO/AN", "INFO/AF",
             "INFO/ALGORITHMS", "INFO/SOFT_FILTERS", "INFO/gnomad_sv_AF",
             "INFO/gnomad_v4.1_sv_AF", "INFO/PREDICTED_LOF", "INFO/PREDICTED_INTERGENIC"]:
    check(f"keeps {want}", want in tags)
for junk in ["INFO/PRPOS", "INFO/PREND", "INFO/SVDB_JUNK"]:
    check(f"drops {junk}", junk not in tags)
check("absent CHR2 not listed", "INFO/CHR2" not in tags,
      "bcftools errors on tags the header lacks")

print("FORMAT keep-list contents")
check("keeps FORMAT/GT (keep-mode)", "^FORMAT/GT" in tags)
for want in ["FORMAT/DHFFC", "FORMAT/DHBFC", "FORMAT/DHFC", "FORMAT/DHBZ"]:
    check(f"keeps {want}", want in tags)
for junk in ["FORMAT/PL", "FORMAT/PR", "FORMAT/SR"]:
    check(f"drops {junk}", junk not in tags and "^" + junk not in tags)

print("no FORMAT declarations means no FORMAT clause, not a bare ^FORMAT")
r = run_keep("""##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
""")
check("awk exits 0", r.returncode == 0, r.stderr[:200])
check("no FORMAT entry", "FORMAT" not in r.stdout, r.stdout[:120])

print("empty keep-list is a hard failure, not a bare ^")
r = run_keep("""##fileformat=VCFv4.2
##INFO=<ID=PRPOS,Number=.,Type=String,Description="junk only">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
""")
check("nonzero exit", r.returncode != 0)
check("no stdout", r.stdout.strip() == "", r.stdout[:120])

print("bcftools round-trip")
if shutil.which("bcftools"):
    with tempfile.TemporaryDirectory() as td:
        vcf = Path(td) / "in.vcf"
        vcf.write_text(HEADER + "\n" + RECORD + "\n")
        keep = subprocess.run(["awk", "-v", "GNOMAD_POP=gnomad_v4.1", "-f", str(KEEP),
                               str(vcf)], capture_output=True, text=True).stdout.strip()
        out_vcf = Path(td) / "out.vcf"
        r = subprocess.run(["bcftools", "annotate", "-x", keep, "-o", str(out_vcf),
                            str(vcf)], capture_output=True, text=True)
        check("bcftools annotate exits 0", r.returncode == 0, r.stderr[:300])
        stripped = out_vcf.read_text()
        check("PRPOS gone from records", "PRPOS=" not in stripped)
        check("PRPOS header line gone", "ID=PRPOS" not in stripped)
        check("SVDB_JUNK gone", "SVDB_JUNK" not in stripped)
        check("SVTYPE survives on record", "SVTYPE=DEL" in stripped)
        check("PREDICTED_LOF survives", "PREDICTED_LOF=GENE1" in stripped)
        check("gnomad_sv_AF survives", "gnomad_sv_AF=0.001" in stripped)
        data = stripped.rstrip().split("\n")[-1].split("\t")
        check("FORMAT stripped to GT + duphold", data[8] == "GT:DHFFC:DHBFC:DHFC:DHBZ",
              data[8])
        check("SR header line gone", "ID=SR" not in stripped)
        check("comma-valued SR gone from sample", "98,0" not in stripped)
        check("GT and duphold values survive", data[9] == "0/1:0.61:0.58:0.63:-1.2",
              data[9])
else:
    print("  skip  bcftools not on PATH -- run inside the test env for the round-trip")

print()
if failures:
    print(f"{len(failures)} failure(s): {', '.join(failures)}")
    sys.exit(1)
print("all checks passed")
