#!/usr/bin/env python3
"""Regression tests for the AnnotSV TSV filter tiers.

Covers nf/annotate_svs/assets/filter_annotsv_tsv.awk directly -- no Nextflow, no container,
and no AnnotSV, which needs a multi-gigabyte annotation bundle.

This is the one script in the pipeline that REMOVES rows, so a wrong column index or an
inverted comparison does not error -- it silently decides which variants a reviewer ever
sees. The column-order tests are the ones to read first: AnnotSV's column set moves between
releases, and a filter reading the wrong column is worse than no filter at all.

Read the row-type tests second. AnnotSV runs in `-annotationMode both`, so the unfiltered TSV
holds full rows (per SV) and split rows (per SV x gene). Both tiers keep split rows only,
because every gene-level column a criterion reads is blank on a full row -- a filter that
accepted full rows would duplicate every variant and evaluate coding-ness against blanks.

    python tests/test_tsv_filters.py

Needs `awk` only.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

AWK = (Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
       / "filter_annotsv_tsv.awk")

COLUMNS = ["AnnotSV_ID", "ID", "FILTER", "Annotation_mode", "Location",
           "Overlapped_CDS_percent", "ACMG_class", "GnomAD_pLI", "LOEUF_bin",
           "B_gain_AFmax", "B_loss_AFmax", "OMIM_morbid", "AnnotSV_ranking_score"]

# A row that passes tier 1 on every criterion. Each test perturbs one field from here, so a
# failure names the single criterion that moved.
PASSING = {
    "AnnotSV_ID": "sv_1", "ID": "rec_1", "FILTER": "PASS", "Annotation_mode": "split",
    "Location": "CDS12-CDS20", "Overlapped_CDS_percent": "43", "ACMG_class": "4",
    "GnomAD_pLI": "0.97", "LOEUF_bin": "1", "B_gain_AFmax": "", "B_loss_AFmax": "",
    "OMIM_morbid": "yes", "AnnotSV_ranking_score": "0.95",
}

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def run(rows, tier=1, columns=COLUMNS, **kw):
    """Filter `rows` (dicts of column -> value) and return (result, kept ids)."""
    lines = ["\t".join(columns)]
    for r in rows:
        lines.append("\t".join(str(r.get(col, "")) for col in columns))
    args = ["awk", "-v", f"TIER={tier}"]
    for k, v in kw.items():
        args += ["-v", f"{k}={v}"]
    args += ["-f", str(AWK)]
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "annotsv.tsv"
        src.write_text("\n".join(lines) + "\n")
        r = subprocess.run(args + [str(src)], capture_output=True, text=True)
    kept = []
    id_at = columns.index("ID")
    for line in r.stdout.splitlines():
        if line.startswith("##") or line.startswith(columns[0] + "\t"):
            continue
        kept.append(line.split("\t")[id_at])
    return r, kept


def row(**overrides):
    r = dict(PASSING)
    r.update(overrides)
    return r


def main():
    if not AWK.is_file():
        sys.exit(f"missing {AWK}")

    print("the filtered file says what it filtered on")
    r, kept = run([row()])
    head = r.stdout.splitlines()[0]
    check("first line is a comment", head.startswith("##"), head[:80])
    check("names the tier", "tier1" in head, head[:80])
    check("carries the active thresholds, not just their names",
          "0.9" in head and "0.01" in head, head[:120])
    check("the AnnotSV header follows it verbatim",
          r.stdout.splitlines()[1] == "\t".join(COLUMNS))
    check("a passing row survives", kept == ["rec_1"], r.stdout)

    r, _ = run([row()], RARE_AF="0.05", PLI_MIN="0.5")
    check("overridden thresholds reach the header line",
          "0.05" in r.stdout.splitlines()[0] and "0.5" in r.stdout.splitlines()[0],
          "six months from now the file has to answer this itself")

    print("tier 1 -- every criterion actually bites")
    for name, over in [
        ("a population tag", {"FILTER": "COMMON_GNOMAD"}),
        ("a caller-support tag", {"FILTER": "NO_CALLER_SUPPORT"}),
        ("a depth tag", {"FILTER": "DEPTH_UNSUPPORTED"}),
        ("a full row", {"Annotation_mode": "full"}),
        ("non-coding", {"Location": "intron3-intron4", "Overlapped_CDS_percent": "0"}),
        ("ACMG 3", {"ACMG_class": "3"}),
        ("an unconstrained gene", {"GnomAD_pLI": "0.1", "LOEUF_bin": "8"}),
        ("a common benign gain", {"B_gain_AFmax": "0.4"}),
        ("a common benign loss", {"B_loss_AFmax": "0.4"}),
    ]:
        _, kept = run([row(**over)])
        check(f"drops {name}", kept == [], f"kept {kept}")

    _, kept = run([row(ACMG_class="5")])
    check("keeps ACMG 5", kept == ["rec_1"])
    _, kept = run([row(GnomAD_pLI="0.1")])
    check("a low LOEUF bin alone is enough", kept == ["rec_1"],
          "constraint is pLI OR LOEUF, not both")
    _, kept = run([row(LOEUF_bin="9")])
    check("a high pLI alone is enough", kept == ["rec_1"])

    print("blank benign AF is the rare case, so it passes")
    _, kept = run([row(B_gain_AFmax="", B_loss_AFmax="")])
    check("no overlapping benign SV is not 'common'", kept == ["rec_1"],
          "treating blank as failing would empty the file on exactly the interesting rows")

    print("split rows only -- the whole reason AnnotSV runs in `both` mode")
    # A full row is the same SV without its gene-level annotation: Location,
    # Overlapped_CDS_percent, GnomAD_pLI and OMIM_morbid are all blank on it. Keeping it would
    # duplicate every variant AND evaluate the coding criterion against nothing.
    full = row(ID="full_row", Annotation_mode="full", Location="",
               Overlapped_CDS_percent="", GnomAD_pLI="", LOEUF_bin="", OMIM_morbid="")
    split = row(ID="split_row")
    for tier in (1, 2):
        _, kept = run([full, split], tier=tier)
        check(f"tier {tier} keeps the split row and drops its full parent",
              kept == ["split_row"], f"kept {kept}")

    _, kept = run([row(Location="", Overlapped_CDS_percent="")])
    check("a split row with no CDS evidence is non-coding, not 'unknown'", kept == [],
          "on a split row the columns are populated, so blank is a real negative")
    _, kept = run([row(Location="", Overlapped_CDS_percent="0")])
    check("an explicit zero CDS overlap drops", kept == [])
    _, kept = run([row(Location="intron3-intron4", Overlapped_CDS_percent="12")])
    check("CDS percent alone is enough", kept == ["rec_1"],
          "Location naming introns does not override a measured CDS overlap")

    print("multi-gene rows are scored on their worst gene")
    _, kept = run([row(GnomAD_pLI="0.01;0.97;0.2", LOEUF_bin="9;8")])
    check("the highest pLI in the list counts", kept == ["rec_1"])
    _, kept = run([row(GnomAD_pLI="0.01;0.02", LOEUF_bin="9;1")])
    check("the lowest LOEUF bin in the list counts", kept == ["rec_1"])
    _, kept = run([row(B_gain_AFmax="0.0001;0.4")])
    check("the highest benign AF in the list counts", kept == [],
          "taking the first value would let a common overlap through")

    print("ACMG_class is read leniently")
    _, kept = run([row(ACMG_class="full=4")])
    check("a 'full=4' cell is class 4", kept == ["rec_1"])
    _, kept = run([row(ACMG_class="NA")])
    check("an NA cell is not a pass", kept == [])

    print("tier 2 -- the working list")
    common = row(ID="common", FILTER="COMMON_GNOMAD")
    internal = row(ID="internal", FILTER="COMMON_INTERNAL")
    nocaller = row(ID="nocaller", FILTER="NO_CALLER_SUPPORT")
    depth = row(ID="depth", FILTER="DEPTH_UNSUPPORTED;NO_CALLER_SUPPORT")
    vus = row(ID="vus", ACMG_class="3", GnomAD_pLI="0.01", LOEUF_bin="9")
    _, kept = run([common, internal, nocaller, depth, vus], tier=2)
    check("drops the population-common rows",
          "common" not in kept and "internal" not in kept, f"kept {kept}")
    check("keeps NO_CALLER_SUPPORT", "nocaller" in kept,
          "a caller missing an event is weak evidence of absence")
    check("keeps DEPTH_UNSUPPORTED", "depth" in kept)
    check("keeps a VUS in an unconstrained gene", "vus" in kept,
          "tier 2 does not apply the constraint criterion")

    r, kept = run([nocaller], tier=2)
    check("the tag is still visible in the output", "NO_CALLER_SUPPORT" in r.stdout,
          "filtering on a tag must not hide it from the reviewer")

    print("tier 2 keeps regulatory hits in disease genes")
    promoter = row(ID="promoter", Location="txStart-txEnd",
                   Overlapped_CDS_percent="0", OMIM_morbid="yes", ACMG_class="3")
    _, kept = run([promoter], tier=2)
    check("non-coding but OMIM morbid survives", kept == ["promoter"],
          "coding OR disease gene is the whole point of the second tier")

    noncoding = row(ID="noncoding", Location="txStart-txEnd", Overlapped_CDS_percent="0",
                    OMIM_morbid="", ACMG_class="3", AnnotSV_ranking_score="0.1")
    _, kept = run([noncoding], tier=2)
    check("non-coding outside a disease gene drops", kept == [], f"kept {kept}")

    print("tier 2's ranking score rescues rows AnnotSV could not class")
    unclassed = row(ID="unclassed", ACMG_class="NA", AnnotSV_ranking_score="0.95")
    _, kept = run([unclassed], tier=2)
    check("a high score is enough on its own", kept == ["unclassed"])
    _, kept = run([row(ID="low", ACMG_class="NA", AnnotSV_ranking_score="0.1")], tier=2)
    check("a low score is not", kept == [])
    _, kept = run([unclassed], tier=2, RANK_MIN="0.99")
    check("the threshold is a knob", kept == [])

    print("columns are found by name, not position")
    shuffled = list(reversed(COLUMNS))
    _, kept = run([row()], columns=shuffled)
    check("reversed column order still passes the same row", kept == ["rec_1"],
          "AnnotSV's column set moves between releases")
    extra = ["Some_New_Column"] + COLUMNS + ["Another_One"]
    _, kept = run([row()], columns=extra)
    check("added columns do not shift the reads", kept == ["rec_1"])
    _, kept = run([row(FILTER="COMMON_GNOMAD")], columns=shuffled)
    check("and a failing row still fails when shuffled", kept == [],
          "a filter that reads the wrong column is worse than no filter")

    print("a renamed column is reported, not guessed at")
    renamed = [c if c != "GnomAD_pLI" else "gnomAD_pLI_v5" for c in COLUMNS]
    r, _ = run([row()], columns=renamed)
    head = r.stdout.splitlines()[0]
    check("the header line names it",
          "MISSING COLUMNS" in head and "GnomAD_pLI" in head, head[:160])

    print("empty input")
    r, kept = run([])
    check("exits 0", r.returncode == 0, r.stderr.strip()[:120])
    check("still writes the provenance line and the header",
          len(r.stdout.splitlines()) == 2, r.stdout)

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("tsv filtering holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
