#!/usr/bin/env python3
"""Regression tests for the somalier -> ploidy table.

Covers nf/shared/assets/check_somalier.awk directly -- no Nextflow, no container, no
somalier.

    python3 tests/test_check_somalier.py

WHY THIS EXISTS. Everything downstream of this file trusts it: the depth filter picks a
chrX/chrY threshold from x_copies/y_copies, and the Talos tail's sex-stratified AF reads the
PED this writes. Both failures are silent -- a sample called female when it is male gets its
chrX duplications judged against the diploid threshold and quietly loses them, and a wrong
column-5 value shifts an AF a recessive filter then acts on.

The PED is a cross-check, never the source: where the two disagree the DATA wins
and the disagreement is made loud. That direction is asserted here, because reversing it
would still look like a working pipeline.

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

AWK = (Path(__file__).parent.parent / "nf" / "shared" / "assets"
       / "check_somalier.awk")

SAMPLES_HEADER = "#family_id\tsample_id\tdepth_mean\tX_depth_mean\tY_depth_mean"

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run(samples_rows, ped_rows, pairs_rows=None, **kw):
    """Run the awk over synthetic tables; return (ploidy rows, inferred ped, stderr, rc).

    ploidy rows come back as a dict keyed on sample so a test names the sample it means
    rather than a row index.
    """
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "samples.tsv").write_text(
            SAMPLES_HEADER + "\n" + "".join(r + "\n" for r in samples_rows))
        (tmp / "cohort.ped").write_text("".join(r + "\n" for r in ped_rows))
        args = ["awk",
                "-v", f"SAMPLES={tmp}/samples.tsv",
                "-v", f"PED={tmp}/cohort.ped",
                "-v", f"OUT_PLOIDY={tmp}/ploidy.tsv",
                "-v", f"OUT_PED={tmp}/cohort.inferred.ped"]
        if pairs_rows is not None:
            (tmp / "pairs.tsv").write_text("".join(r + "\n" for r in pairs_rows))
            args += ["-v", f"PAIRS={tmp}/pairs.tsv"]
        for k, v in kw.items():
            args += ["-v", f"{k}={v}"]
        args += ["-f", str(AWK)]
        r = subprocess.run(args, stdin=subprocess.DEVNULL, capture_output=True, text=True)

        ploidy, inferred = {}, []
        if (tmp / "ploidy.tsv").exists():
            for line in (tmp / "ploidy.tsv").read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split("\t")
                ploidy[f[0]] = f[1:]
        if (tmp / "cohort.inferred.ped").exists():
            inferred = (tmp / "cohort.inferred.ped").read_text().splitlines()
    return ploidy, inferred, r.stderr, r.returncode


# A male: X at half autosomal depth, Y present. A female: X at full depth, Y at noise.
MALE = "F1\tS1\t30\t15\t12"
FEMALE = "F1\tS2\t30\t30\t0.6"
PED_MALE = "F1\tS1\t0\t0\t1\t2"
PED_FEMALE = "F1\tS2\t0\t0\t2\t1"


def test_karyotypes():
    print("karyotype calls")
    ploidy, _, err, rc = run([MALE, FEMALE], [PED_MALE, PED_FEMALE])
    check("exit 0 on a clean cohort", rc == 0, err.strip()[:200])
    check("male is 1 X, 1 Y", ploidy.get("S1", [])[:2] == ["1", "1"], str(ploidy))
    check("female is 2 X, 0 Y", ploidy.get("S2", [])[:2] == ["2", "0"], str(ploidy))
    check("male agrees with its PED row", ploidy["S1"][3] == "AGREES", str(ploidy["S1"]))
    check("female agrees with its PED row", ploidy["S2"][3] == "AGREES", str(ploidy["S2"]))

    # X ratio 0.72 sits in the gap between the hemizygous and diploid cut points. Rounding it
    # to the nearer mode is exactly the confident-wrong answer this column exists to avoid.
    ploidy, _, err, rc = run(["F1\tS3\t30\t21.6\t6"], ["F1\tS3\t0\t0\t1\t2"])
    check("ambiguous X ratio is undetermined", ploidy.get("S3", [])[0] == ".", str(ploidy))
    check("undetermined karyotype reports UNKNOWN",
          ploidy.get("S3", [])[3] == "UNKNOWN", str(ploidy))
    check("undetermined karyotype is noted on stderr", "no determined karyotype" in err, err)

    # XXY: two X, one Y. Per-chromosome copies are still usable even though "sex" is not.
    ploidy, _, err, rc = run(["F1\tS4\t30\t30\t12"], ["F1\tS4\t0\t0\t1\t2"])
    check("XXY keeps usable per-chromosome copies",
          ploidy.get("S4", [])[:2] == ["2", "1"], str(ploidy))
    check("XXY has no inferred sex", ploidy.get("S4", [])[3] == "UNKNOWN", str(ploidy))


def test_ped_disagreement():
    print("PED cross-check")
    # The PED says female, the alignments say male.
    ploidy, inferred, err, rc = run([MALE], ["F1\tS1\t0\t0\t2\t2"])
    check("disagreement is flagged in the table",
          ploidy.get("S1", [])[3] == "DISAGREES", str(ploidy))
    check("disagreement warns on stderr", "WARNING sample S1" in err, err)
    check("the run continues", rc == 0, err.strip()[:200])
    check("the DATA wins in the inferred PED",
          inferred and inferred[0].split("\t")[4] == "1", str(inferred))

    ploidy, inferred, err, _ = run([MALE], ["F1\tS1\t0\t0\t0\t2"])
    check("PED sex 0 reads as missing, not as a disagreement",
          ploidy.get("S1", [])[3] == "PED_MISSING", str(ploidy))
    check("a missing PED sex is filled in from the data",
          inferred[0].split("\t")[4] == "1", str(inferred))


def test_inferred_ped():
    print("inferred PED")
    ped = ["#family sample father mother sex phenotype",
           "F1\tS1\tDAD\tMUM\t2\t2",
           "F1\tS3\t0\t0\t1\t1"]
    # S1 is male in the data; S3 is ambiguous, so its PED row must survive untouched.
    _, inferred, err, rc = run([MALE, "F1\tS3\t30\t21.6\t6"], ped)
    check("comment lines survive", inferred[0].startswith("#"), str(inferred))
    row = inferred[1].split("\t")
    check("column 5 is replaced with the inferred sex", row[4] == "1", str(row))
    check("family and parent columns are untouched",
          row[0] == "F1" and row[2] == "DAD" and row[3] == "MUM", str(row))
    check("phenotype is untouched", row[5] == "2", str(row))
    check("an undetermined sample keeps its PED sex",
          inferred[2].split("\t")[4] == "1", str(inferred))
    # A PED naming a sample somalier never saw -- one site-wide PED reused across cohorts is
    # normal -- must pass through rather than being dropped or blanked.
    _, inferred, _, rc = run([MALE], [PED_MALE, "F9\tSX\t0\t0\t2\t1"])
    check("a PED row for an unseen sample survives",
          len(inferred) == 2 and inferred[1].split("\t")[1] == "SX", str(inferred))


def test_relatedness():
    print("relatedness cross-check")
    pairs_header = "#sample_a\tsample_b\trelatedness\texpected_relatedness"
    ploidy, _, err, rc = run(
        [MALE, FEMALE], [PED_MALE, PED_FEMALE],
        [pairs_header, "S1\tS2\t0.02\t0.5"])
    check("a declared relation that does not measure warns",
          "declared related" in err, err)
    check("the run continues anyway", rc == 0, err.strip()[:200])

    _, _, err, _ = run([MALE, FEMALE], [PED_MALE, PED_FEMALE],
                       [pairs_header, "S1\tS2\t0.48\t0.0"])
    check("an unexpected relation warns", "declared unrelated" in err, err)

    _, _, err, _ = run([MALE, FEMALE], [PED_MALE, PED_FEMALE],
                       [pairs_header, "S1\tS2\t0.50\t0.5"])
    check("a matching pair is silent", "WARNING" not in err, err)

    # -1 is what somalier writes when the PED implies no expectation for a pair.
    _, _, err, _ = run([MALE, FEMALE], [PED_MALE, PED_FEMALE],
                       [pairs_header, "S1\tS2\t0.60\t-1"])
    check("no expectation means no opinion", "WARNING" not in err, err)


def test_hard_failures():
    print("failures that must be loud")
    # A sites file built against the wrong reference genotypes nothing and would otherwise
    # produce a complete, plausible, entirely wrong ploidy table.
    _, _, err, rc = run(["F1\tS1\t0\t0\t0"], [PED_MALE])
    check("zero depth is fatal", rc != 0, f"rc={rc}")
    check("zero depth names the likely cause", "--somalier_sites" in err, err)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "samples.tsv").write_text("#family_id\tsample_id\tdepth_mean\nF1\tS1\t30\n")
        (tmp / "cohort.ped").write_text(PED_MALE + "\n")
        r = subprocess.run(
            ["awk", "-v", f"SAMPLES={tmp}/samples.tsv", "-v", f"PED={tmp}/cohort.ped",
             "-v", f"OUT_PLOIDY={tmp}/p.tsv", "-v", f"OUT_PED={tmp}/i.ped",
             "-f", str(AWK)],
            stdin=subprocess.DEVNULL, capture_output=True, text=True)
    check("a renamed/absent column is fatal", r.returncode != 0, f"rc={r.returncode}")
    check("the missing column is named", "X_depth_mean" in r.stderr, r.stderr)

    _, _, err, rc = run([MALE], ["F1\tS1\t0\t0\t1"])
    check("a short PED line is fatal", rc != 0, f"rc={rc}")


def main():
    for t in (test_karyotypes, test_ped_disagreement, test_inferred_ped,
              test_relatedness, test_hard_failures):
        t()
    print()
    if failures:
        print(f"{len(failures)} failure(s): {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
