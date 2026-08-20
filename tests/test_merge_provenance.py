#!/usr/bin/env python3
"""Regression tests for the between-axes provenance promotion.

Covers nf/annotate_svs/assets/promote_caller_support.awk directly -- no Nextflow, no
container, so it gates every change to the logic that keeps the two merge stages from
colliding.

Why this exists: the design spike established by experiment that axis B carries
axis A's INFO keys through verbatim and then appends its own, leaving two svdb_origin, two
FOUNDBY, two SUPP_VEC, two set and two VARID on a twice-merged record -- with bcftools
returning the stale axis-A value. Nothing errors. The cohort VCF just quietly reports the
wrong support counts. This awk strips the collision before axis B can create it.

    python tests/test_merge_provenance.py

Needs `awk` only.
"""

import re
import subprocess
import sys
import tempfile
from pathlib import Path

AWK = (Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
       / "promote_caller_support.awk")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1"""

# What SVDB actually emits for a two-caller axis-A merge, taken from the design spike output.
MERGED = ("SVTYPE=DEL;END=2000;ALGORITHMS=manta;CALLER_SUPP=manta;NCALLER=1;"
          "VARID=delly_1:delly;set=Intersection;FOUNDBY=2;"
          "manta_INFO=manta_1|SVTYPE:DEL;delly_INFO=delly_1|SVTYPE:DEL;"
          "svdb_origin=manta|delly;SUPP_VEC=11")

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def promote(info_strings):
    return _run(info_strings, [])


def _run(info_strings, extra_args):
    vcf = HEADER + "\n" + "\n".join(
        f"chr1\t{1000 + i}\trec_{i}\tN\t<DEL>\t99\tPASS\t{info}\tGT\t0/1"
        for i, info in enumerate(info_strings)) + "\n"
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.vcf"
        src.write_text(vcf)
        # Named twice: pass 1 collects the contigs the records use, so the header can
        # declare any SVDB dropped. Once prints nothing -- see the awk's usage note.
        res = subprocess.run(["awk"] + extra_args + ["-f", str(AWK), str(src), str(src)],
                             capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"awk failed: {res.stderr.strip()[:300]}")
    lines = res.stdout.splitlines()
    return ([l for l in lines if l.startswith("#")],
            [l.split("\t")[7] for l in lines if not l.startswith("#")])


def promote_cohort(info_strings):
    """Same as promote() but in axis-B mode."""
    return _run(info_strings, ["-v", "MODE=cohort"])


def keys(info):
    return [kv.split("=")[0] for kv in info.split(";") if "=" in kv]


def val(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("axis-A provenance is promoted to typed keys")
    _, infos = promote([MERGED])
    info = infos[0]
    check("svdb_origin becomes CALLER_SUPP, pipe-separated to comma",
          val(info, "CALLER_SUPP") == "manta,delly",
          f"got {val(info, 'CALLER_SUPP')!r}")
    check("FOUNDBY becomes NCALLER", val(info, "NCALLER") == "2",
          f"got {val(info, 'NCALLER')!r}")
    check("a stale pre-merge CALLER_SUPP is replaced, not kept alongside",
          keys(info).count("CALLER_SUPP") == 1 and val(info, "CALLER_SUPP") != "manta",
          f"INFO={info}")

    print("the raw merge keys are gone before axis B can duplicate them")
    for key in ("svdb_origin", "FOUNDBY", "SUPP_VEC", "set", "VARID"):
        check(f"{key} is stripped", key not in keys(info),
              "left in INFO -- axis B would append its own and bcftools would "
              "return this stale one")
    check("no INFO key appears twice",
          len(keys(info)) == len(set(keys(info))), f"keys={keys(info)}")

    print("unrelated INFO survives")
    check("SVTYPE survives", val(info, "SVTYPE") == "DEL", f"INFO={info}")
    check("END survives", val(info, "END") == "2000", f"INFO={info}")
    check("per-input <tag>_INFO survives",
          val(info, "manta_INFO") is not None and val(info, "delly_INFO") is not None,
          "these carry each input's original record, worth keeping")

    print("ALGORITHMS is re-derived from the support, not left at the priority caller")
    check("ALGORITHMS names every contributing caller",
          val(info, "ALGORITHMS") == "manta,delly",
          f"got {val(info, 'ALGORITHMS')!r} -- SVDB keeps only the priority record's "
          "value, so unfixed this reports one caller where two agreed")
    check("ALGORITHMS appears once", keys(info).count("ALGORITHMS") == 1,
          f"keys={keys(info)}")

    print("records that never went through a merge are left alone")
    solo = "SVTYPE=DEL;END=2000;ALGORITHMS=tiddit;CALLER_SUPP=tiddit;NCALLER=1"
    _, infos = promote([solo])
    info = infos[0]
    check("CALLER_SUPP from standardization is preserved",
          val(info, "CALLER_SUPP") == "tiddit", f"INFO={info}")
    check("NCALLER from standardization is preserved",
          val(info, "NCALLER") == "1", f"INFO={info}")
    check("ALGORITHMS is untouched when there is no merge provenance",
          val(info, "ALGORITHMS") == "tiddit", f"INFO={info}")

    print("edge cases")
    _, infos = promote(["svdb_origin=manta|delly|wham;FOUNDBY=3"])
    check("three callers promote cleanly",
          val(infos[0], "CALLER_SUPP") == "manta,delly,wham"
          and val(infos[0], "NCALLER") == "3", f"INFO={infos[0]}")

    _, infos = promote(["svdb_origin=manta"])
    check("NCALLER is derived when FOUNDBY is absent",
          val(infos[0], "NCALLER") == "1", f"INFO={infos[0]}")

    _, infos = promote(["SUPP_VEC=11;set=Intersection"])
    check("stripping every key leaves a valid '.' INFO",
          infos[0] == ".", f"got {infos[0]!r}")

    hdr, _ = promote([MERGED])
    for key in ("CALLER_SUPP", "NCALLER"):
        check(f"{key} is declared exactly once in the header",
              sum(1 for l in hdr if l.startswith(f"##INFO=<ID={key}")) == 1,
              "bcftools needs it declared to read it typed, and twice is malformed")

    print("MODE=cohort recovers callers from the per-input blobs, not from sample sets")
    # At axis B svdb_origin names SAMPLE SETS. Promoting it put CALLER_SUPP=COHORT,SAMP1,
    # SAMP3 on the cohort VCF -- sample-set labels masquerading as caller support, which is
    # exactly the conflation to avoid. The callers come from the <tag>_INFO blobs instead.
    # Note SVDB rewrites a nested comma as a colon, hence "tiddit:cnvpytor" below.
    cohort_info = (
        "SVTYPE=DEL;END=2000;ALGORITHMS=dragen;CALLER_SUPP=dragen;NCALLER=1;"
        "SAMP1_INFO=x|SVTYPE:DEL|CALLER_SUPP:tiddit:cnvpytor|NCALLER:2;"
        "COHORT_INFO=y|SVTYPE:DEL|CALLER_SUPP:dragen|NCALLER:1;"
        "SAMP3_INFO=z|SVTYPE:DEL|CALLER_SUPP:tiddit|NCALLER:1;"
        "svdb_origin=COHORT|SAMP1|SAMP3;SUPP_VEC=111;FOUNDBY=3")
    _, infos = promote_cohort([cohort_info])
    info = infos[0]
    supp = (val(info, "CALLER_SUPP") or "").split(",")
    check("CALLER_SUPP is the union of the inputs' callers",
          sorted(supp) == ["cnvpytor", "dragen", "tiddit"],
          f"got {val(info, 'CALLER_SUPP')!r} -- sample-set labels here means axis B's "
          "svdb_origin was promoted by mistake")
    check("each caller appears once despite two inputs reporting tiddit",
          supp.count("tiddit") == 1, f"got {val(info, 'CALLER_SUPP')!r}")
    check("NCALLER counts distinct callers", val(info, "NCALLER") == "3",
          f"got {val(info, 'NCALLER')!r}")
    check("ALGORITHMS mirrors the union", val(info, "ALGORITHMS") == val(info, "CALLER_SUPP"),
          f"ALG={val(info, 'ALGORITHMS')!r} SUPP={val(info, 'CALLER_SUPP')!r}")
    check("axis B's own svdb_origin is kept as sample-set provenance",
          val(info, "svdb_origin") == "COHORT|SAMP1|SAMP3",
          "it names which sample sets contributed -- worth keeping, and there is no "
          "later merge for it to collide with")
    check("axis B's SUPP_VEC is kept", val(info, "SUPP_VEC") == "111", f"INFO={info}")
    check("no INFO key appears twice at cohort level",
          len(keys(info)) == len(set(keys(info))), f"keys={keys(info)}")

    # A cohort record from a single input still has to report its caller, not its label.
    solo_cohort = ("SVTYPE=DEL;CALLER_SUPP=SAMP1;NCALLER=1;"
                   "SAMP1_INFO=x|SVTYPE:DEL|CALLER_SUPP:tiddit|NCALLER:1;"
                   "svdb_origin=SAMP1;SUPP_VEC=010")
    _, infos = promote_cohort([solo_cohort])
    check("a single-input cohort record reports its caller, not its sample set",
          val(infos[0], "CALLER_SUPP") == "tiddit" and val(infos[0], "NCALLER") == "1",
          f"INFO={infos[0]}")

    print("a semicolon inside a per-input blob does not split the INFO field")
    # SVDB copies a source record's sample and FILTER columns into <tag>_SAMPLE/_FILTERS
    # verbatim and does not escape every semicolon it copies. INFO is semicolon-delimited,
    # so one such value ends the field early and the remainder is not KEY=VALUE -- bcftools
    # refuses the record ("not defined in the header") when promote_support's output is read
    # back by bcftools sort. That is what killed the first real cohort run at axis B.
    stray = ("SVTYPE=DEL;END=2000;"
             "SAMP1_SAMPLE=x|SAMP1|GT:0/1;PR:5,3|SAMP2|GT:1/1;"
             "SAMP1_INFO=x|SVTYPE:DEL|CALLER_SUPP:manta|NCALLER:1;"
             "svdb_origin=SAMP1;SUPP_VEC=1")
    _, infos = promote_cohort([stray])
    info = infos[0]
    check("every INFO field is KEY=VALUE or a bare flag",
          all(re.match(r"^[A-Za-z_][A-Za-z0-9_.]*(=|$)", kv) for kv in info.split(";")),
          f"INFO={info} -- a fragment here is what bcftools cannot resolve")
    check("the orphaned fragment is reattached to the blob it came from",
          val(info, "SAMP1_SAMPLE") == "x|SAMP1|GT:0/1|PR:5,3|SAMP2|GT:1/1",
          f"got {val(info, 'SAMP1_SAMPLE')!r}")
    check("the blob after the break is still read for caller support",
          val(info, "CALLER_SUPP") == "manta", f"INFO={info}")

    _, infos = promote(["SVTYPE=DEL;IMPRECISE;manta_INFO=x|SVTYPE:DEL;svdb_origin=manta"])
    check("a genuine INFO flag is not swallowed by the field before it",
          "IMPRECISE" in infos[0].split(";"), f"INFO={infos[0]}")

    _, infos = promote(["SVTYPE=DEL;manta_INFO=x|FILTER:a;b;svdb_origin=manta|delly"])
    check("axis A gets the same repair",
          all(re.match(r"^[A-Za-z_][A-Za-z0-9_.]*(=|$)", kv)
              for kv in infos[0].split(";"))
          and val(infos[0], "manta_INFO") == "x|FILTER:a|b",
          f"INFO={infos[0]}")

    print("contigs the header never declared are backfilled")
    # SVDB copies one input's header into the merged VCF, so a contig only the other input
    # called on goes undeclared. bcftools sort then dies -- it flushes a temp chunk written
    # against the header it was handed, htslib auto-adds the contig mid-stream, and the two
    # disagree ("Error encountered while parsing the input at chr19:1", seen on real data
    # where cnvnator called chr19 and the kept header did not name it).
    vcf = (HEADER + "\n"
           + "chr1\t1000\tr0\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL\tGT\t0/1\n"
           + "chr19\t1\tr1\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL\tGT\t1/1\n")
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.vcf"
        src.write_text(vcf)
        out = subprocess.run(["awk", "-f", str(AWK), str(src), str(src)],
                             capture_output=True, text=True).stdout
    hdr = [l for l in out.splitlines() if l.startswith("##contig=")]
    check("the undeclared contig gains a header line",
          any(l.startswith("##contig=<ID=chr19") for l in hdr), f"contigs={hdr}")
    check("the already-declared contig is not declared twice",
          sum(1 for l in hdr if l.startswith("##contig=<ID=chr1,")
              or l == "##contig=<ID=chr1>") == 1, f"contigs={hdr}")
    check("both records survive the extra pass",
          len([l for l in out.splitlines() if not l.startswith("#")]) == 2,
          f"got {len([l for l in out.splitlines() if not l.startswith('#')])} records")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("provenance promotion holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
