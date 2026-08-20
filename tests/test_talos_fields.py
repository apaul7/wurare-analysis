#!/usr/bin/env python3
"""Regression tests for the Talos schema transform and field report.

Covers nf/annotate_svs/assets/talos_schema.awk and check_talos_fields.awk directly -- no
Nextflow, no container, no Talos.

Talos's failure mode is silence, which is what makes these worth having. Two documented
examples from the design: a variant with no PREDICTED_LOF entry is dropped
entirely rather than reported, and an absent ALGORITHMS is DEFAULTED to ['gCNV'] rather than
rejected. A callset can pass through Talos and come out mislabelled or half empty with
nothing anywhere saying so.

    python tests/test_talos_fields.py

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

ASSETS = Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
SCHEMA = ASSETS / "talos_schema.awk"
CHECK = ASSETS / "check_talos_fields.awk"

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Len">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=AC,Number=1,Type=Integer,Description="AC">
##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">
##INFO=<ID=AF,Number=1,Type=Float,Description="AF">
##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Algs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3"""

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run_awk(awk, vcf_text, pop="gnomad_v4.1"):
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.vcf"
        src.write_text(vcf_text)
        r = subprocess.run(["awk", "-v", f"GNOMAD_POP={pop}", "-f", str(awk), str(src)],
                           capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"awk failed: {r.stderr.strip()[:300]}")
    return r.stdout


def build(info, gts=("0/1", "1/1", "0/0"), extra_info_headers=(), filt="PASS"):
    head = HEADER
    if extra_info_headers:
        lines = head.splitlines()
        decls = [f'##INFO=<ID={i},Number=1,Type=String,Description="x">'
                 for i in extra_info_headers]
        head = "\n".join(lines[:-1] + decls + lines[-1:])
    return (head + "\n"
            + f"chr1\t1000\trec_1\tN\t<DEL>\t99\t{filt}\t{info}\tGT\t" + "\t".join(gts) + "\n")


def record(vcf_text):
    return [l for l in vcf_text.splitlines() if not l.startswith("#")][0].split("\t")


def declared_ids(vcf_text):
    return [l.split("##INFO=<ID=")[1].split(",")[0]
            for l in vcf_text.splitlines() if l.startswith("##INFO=<ID=")]


def schema(info, gts=("0/1", "1/1", "0/0"), pop="gnomad_v4.1"):
    out = run_awk(SCHEMA, build(info, gts), pop)
    rec = [l for l in out.splitlines() if not l.startswith("#")][0]
    return rec.split("\t")[7], out


def val(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def main():
    for a in (SCHEMA, CHECK):
        if not a.is_file():
            sys.exit(f"missing {a}")

    print("carrier counts are sample counts, not allele counts")
    # bcftools +fill-tags -t AC_Hom counts ALLELES: it is twice the homozygote count. Talos
    # wants samples, and a factor-of-two error in a carrier count reads as plausible.
    info, _ = schema("SVTYPE=DEL;AC=3;AN=6;AF=0.5", ("0/1", "1/1", "0/0"))
    check("N_HET counts heterozygous samples", val(info, "N_HET") == "1", f"INFO={info}")
    check("N_HOMALT counts hom-alt samples, not alleles",
          val(info, "N_HOMALT") == "1",
          f"got {val(info, 'N_HOMALT')} -- 2 would mean alleles were counted")
    info, _ = schema("SVTYPE=DEL", ("0/1", "0|1", "1/1"))
    check("phased genotypes are counted too",
          val(info, "N_HET") == "2" and val(info, "N_HOMALT") == "1", f"INFO={info}")
    info, _ = schema("SVTYPE=DEL", ("./.", "0/0", "./."))
    check("no-calls and refs are counted as neither",
          val(info, "N_HET") == "0" and val(info, "N_HOMALT") == "0", f"INFO={info}")

    print("population AF is exposed under the name Talos looks for")
    info, _ = schema("SVTYPE=DEL;gnomad_sv_AF=0.021")
    check("gnomad_sv_AF becomes {GNOMAD_POP}_sv_AF",
          val(info, "gnomad_v4.1_sv_AF") == "0.021", f"INFO={info}")
    check("the pipeline's own name survives alongside it",
          val(info, "gnomad_sv_AF") == "0.021",
          "other consumers still read gnomad_sv_AF; a rename would break them")
    info, _ = schema("SVTYPE=DEL;gnomad_sv_AF=0.5", pop="gnomad_v5")
    check("the prefix is configurable", val(info, "gnomad_v5_sv_AF") == "0.5", f"INFO={info}")
    info, _ = schema("SVTYPE=DEL")
    check("a record with no population AF gains no empty field",
          val(info, "gnomad_v4.1_sv_AF") is None, f"INFO={info}")
    info, _ = schema("SVTYPE=DEL;gnomad_sv_AF=.")
    check("a missing-value AF is not copied as '.'",
          val(info, "gnomad_v4.1_sv_AF") is None, f"INFO={info}")

    print("STATUS is stamped only when absent")
    check("absent STATUS is filled",
          val(schema("SVTYPE=DEL")[0], "STATUS") == "PASS")
    check("an existing STATUS is preserved",
          val(schema("SVTYPE=DEL;STATUS=LOWQUAL")[0], "STATUS") == "LOWQUAL",
          "overwriting a caller's own status would discard real information")

    print("existing fields are never clobbered")
    info, _ = schema("SVTYPE=DEL;AC=3;AN=6;AF=0.5;ALGORITHMS=manta,delly;N_HET=99")
    check("AC/AN/AF from Phase 2 survive",
          [val(info, k) for k in ("AC", "AN", "AF")] == ["3", "6", "0.5"], f"INFO={info}")
    check("ALGORITHMS survives -- Talos defaults it to gCNV when absent",
          val(info, "ALGORITHMS") == "manta,delly", f"INFO={info}")
    check("a pre-existing N_HET is not duplicated",
          info.count("N_HET=") == 1, f"INFO={info}")

    print("FILTER tags are carried into INFO, because Talos treats them as deletion")
    # Talos runs mt.filter_rows(is_missing(filters) | filters.length() == 0) BEFORE any
    # category logic, so the soft filter tags erase the record there. Measured: a 2-record VCF
    # tagged COMMON_INTERNAL loaded as 0 rows. The tags are moved, never dropped.
    out = run_awk(SCHEMA, build("SVTYPE=DEL", filt="COMMON_INTERNAL"))
    rec = record(out)
    check("a tagged record is passed through as PASS", rec[6] == "PASS", f"FILTER={rec[6]}")
    check("and the tag survives in SOFT_FILTERS",
          val(rec[7], "SOFT_FILTERS") == "COMMON_INTERNAL", f"INFO={rec[7]}")
    check("SOFT_FILTERS is declared", "SOFT_FILTERS" in declared_ids(out))

    rec = record(run_awk(SCHEMA, build("SVTYPE=DEL",
                                       filt="COMMON_INTERNAL;NO_CALLER_SUPPORT")))
    check("multiple tags join with commas, not semicolons",
          val(rec[7], "SOFT_FILTERS") == "COMMON_INTERNAL,NO_CALLER_SUPPORT",
          f"a semicolon would split INFO into two fields: {rec[7]}")

    rec = record(run_awk(SCHEMA, build("SVTYPE=DEL", filt="PASS")))
    check("an already-passing record gains no SOFT_FILTERS",
          val(rec[7], "SOFT_FILTERS") is None and rec[6] == "PASS", f"INFO={rec[7]}")
    rec = record(run_awk(SCHEMA, build("SVTYPE=DEL", filt=".")))
    check("an unfiltered '.' is left alone, not rewritten to PASS",
          rec[6] == "." and val(rec[7], "SOFT_FILTERS") is None,
          "'.' means no filters were applied and hail already reads it as missing")

    print("every field Talos reads without checking is declared in the header")
    # Talos's rearrange_annotations() tolerates a missing ALGORITHMS, STATUS, CHR2 and END2 and
    # nothing else -- the rest are `mt.info.X` struct accesses, and hail builds that struct from
    # the header. Undeclared means Talos raises on load, which is a different failure from the
    # silent ones above: it never starts rather than quietly under-reporting.
    _, out = schema("SVTYPE=DEL;gnomad_sv_AF=0.01")
    ids = declared_ids(out)
    for f in ("gnomad_v4.1_sv_AF", "gnomad_v4.1_sv_SVID", "N_HET", "N_HOMALT", "STATUS",
              "AF_MALE", "AF_FEMALE"):
        check(f"{f} is declared", f in ids, f"declared: {ids}")
    check("SVID is declared even though this pipeline cannot populate it",
          "gnomad_v4.1_sv_SVID" in ids,
          "svdb query gives a frequency, not an ID -- a null is a degraded report, "
          "an undeclared field is a crash")
    check("the sex-AF pair is declared whole",
          "AF_MALE" in ids and "AF_FEMALE" in ids,
          "Talos branches on AF_MALE then reads AF_FEMALE off the same struct")

    out = run_awk(SCHEMA, build("SVTYPE=DEL", extra_info_headers=("N_HET", "AF_MALE")))
    ids = declared_ids(out)
    check("an already-declared field is not declared twice",
          ids.count("N_HET") == 1 and ids.count("AF_MALE") == 1,
          f"bcftools keeps the first and warns, so a second declaration is silently lost: {ids}")
    check("a caller's own SVTYPE/SVLEN/END declaration is left alone",
          [ids.count(f) for f in ("SVTYPE", "SVLEN", "END")] == [1, 1, 1], f"declared: {ids}")

    # normalize_records.awk derives SVTYPE/SVLEN/END onto records whose caller may never have
    # declared them, and a caller header is the one input here no process in this repo writes.
    bare = "\n".join(l for l in build("SVTYPE=DEL;SVLEN=-100;END=1100").splitlines()
                     if not l.startswith(("##INFO=<ID=SVTYPE", "##INFO=<ID=SVLEN",
                                          "##INFO=<ID=END")))
    ids = declared_ids(run_awk(SCHEMA, bare + "\n"))
    for f in ("SVTYPE", "SVLEN", "END"):
        check(f"an undeclared {f} from a caller VCF is backfilled", f in ids, f"declared: {ids}")
    bare_no_af = "\n".join(l for l in bare.splitlines()
                           if not l.startswith(("##INFO=<ID=AC", "##INFO=<ID=AN",
                                                "##INFO=<ID=AF")))
    ids = declared_ids(run_awk(SCHEMA, bare_no_af + "\n"))
    check("an absent AC/AN/AF is left absent",
          not any(f in ids for f in ("AC", "AN", "AF")),
          "fill_tags writes those unconditionally at Phase 2; backfilling would hide it "
          f"having not run: {ids}")

    print("the field report says what Talos will and will not read")
    full = ("SVTYPE=DEL;SVLEN=-100;END=1100;ALGORITHMS=manta;STATUS=PASS;"
            "PREDICTED_LOF=GENE1;AC=1;AF=0.1;AN=6;N_HET=1;N_HOMALT=0;"
            "gnomad_v4.1_sv_AF=0.001;gnomad_v4.1_sv_SVID=gnomAD-SV_v4.1_DEL_chr1_1")
    out = run_awk(CHECK, build(full))
    check("a complete record reports no MISSING", "MISSING" not in out, out)
    check("PREDICTED_LOF is counted as the hard gate",
          "records_retainable_by_talos\t1" in out, out)

    out = run_awk(CHECK, build("SVTYPE=DEL;AC=1;AF=0.1;AN=6"))
    check("absent PREDICTED_LOF is reported MISSING",
          any(l.startswith("PREDICTED_LOF\t0\tMISSING") for l in out.splitlines()), out)
    check("and the warning is explicit that Talos would retain nothing",
          "Talos would retain nothing" in out,
          "the whole point is that Talos drops these silently")
    check("absent ALGORITHMS is reported MISSING",
          any(l.startswith("ALGORITHMS\t0\tMISSING") for l in out.splitlines()),
          "Talos defaults it to gCNV rather than rejecting -- mislabelled, not caught")

    print("the report separates the crash case from the silent one")
    # An undeclared field is fatal whatever the records hold; a declared-but-empty one is the
    # silent path. Reporting them in one column would merge two different fixes.
    out = run_awk(CHECK, build("SVTYPE=DEL;AC=1;AF=0.1;AN=6"))
    check("an undeclared fatal field reads UNDECLARED, not MISSING",
          any(l.startswith("gnomad_v4.1_sv_SVID\tUNDECLARED") for l in out.splitlines()), out)
    check("and the warning says Talos raises rather than filters badly",
          "Talos raises on load" in out, out)
    check("a declared field reads declared regardless of record coverage",
          any(l.startswith("SVTYPE\tdeclared") for l in out.splitlines()), out)

    check("an absent sex-AF pair is UNDECLARED",
          any(l.startswith("sex-stratified AF pair\tUNDECLARED") for l in out.splitlines()), out)
    out = run_awk(CHECK, build("SVTYPE=DEL", extra_info_headers=("AF_MALE",)))
    check("half a sex-AF pair is still UNDECLARED",
          any(l.startswith("sex-stratified AF pair\tUNDECLARED") for l in out.splitlines()),
          "a single-sex PED declares one of the two, and Talos reads both")
    out = run_awk(CHECK, build("SVTYPE=DEL", extra_info_headers=("MALE_AF", "FEMALE_AF")))
    check("GATK-SV's own MALE_AF/FEMALE_AF naming satisfies the pair",
          any(l.startswith("sex-stratified AF pair\tdeclared as MALE_AF/FEMALE_AF")
              for l in out.splitlines()), out)

    print("partial coverage is distinguished from absence")
    two = (HEADER + "\n"
           + "chr1\t1000\tr1\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL;PREDICTED_LOF=G1\tGT\t0/1\t0/0\t0/0\n"
           + "chr1\t2000\tr2\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL\tGT\t0/1\t0/0\t0/0\n")
    out = run_awk(CHECK, two)
    check("a field on some records reads 'partial', not 'ok' or 'MISSING'",
          any(l.startswith("PREDICTED_LOF\t1\tpartial") for l in out.splitlines()), out)
    check("only the records with the gate are counted retainable",
          "records_retainable_by_talos\t1" in out, out)

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("Talos schema holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
