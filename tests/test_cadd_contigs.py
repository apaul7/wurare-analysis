#!/usr/bin/env python3
"""Regression test for the CADD contig round-trip, chrM in particular.

    python3 tests/test_cadd_contigs.py

CADD v1.6 GRCh38 is Ensembl-named, so run_cadd strips the `chr` prefix on the way in and
merge_cadd puts it back on the way out. The mitochondrion is the one contig whose name is not
just the prefix: Ensembl calls it MT, GATK GRCh38 calls it chrM. A bare `sed 's/^chr//'` sends
`M`, which CADD's MT-keyed lookups miss, and a bare `sed 's/^/chr/'` would turn a returned MT
into chrMT.

Every one of those failures is SILENT. A CADD fragment with no rows for an interval is a
legitimate result -- merge_cadd's `|| true` exists for exactly that -- so a mitochondrion that
was never scored is indistinguishable from one that had no variants, and a chrMT row simply
joins nothing in ANNOVAR. The result is a blank CADD column for chrM with no error anywhere.

Note what is NOT the defect: the old blunt pair round-trips chrM -> M -> chrM quite happily.
The error is the name handed to CADD in the middle. So the fake CADD.sh here scores only the
Ensembl names its real counterpart knows and silently drops everything else -- a fake that
echoed its input back would pass against the broken code, which is what the first version of
this test did.

Ported from an upstream fix that shipped without a test because its tier-1
stubs the process out. Composed here with this repo's conditional chr_prefix, which upstream
does not have -- so the Ensembl-named case is covered too.

Reuses tests/test_skip_params.py's fake-tool environment. Needs `nextflow`; skips without it.
"""

import gzip
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from test_skip_params import build_env, run_pipeline, write_exec

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


# A bcftools that honours -r, so each split VCF really carries its own contig, and that really
# gzips its output, because run_cadd reads it back with zgrep. The shared fake in
# test_skip_params emits chr1 unconditionally, which cannot show a contig being mapped.
def bcftools_fake(chroms):
    listed = " ".join(f"'{c}'" for c in chroms)
    return r"""#!/bin/sh
sub="$1"; shift
out=""; region=""; prev=""; last=""
for a in "$@"; do
  [ "$prev" = "-o" ] && out="$a"
  [ "$prev" = "-r" ] && region="$a"
  prev="$a"; last="$a"
done
hdr='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
case "$sub" in
  query)     printf '%s\n' """ + listed + r""" ;;
  view|norm|concat|sort)
    if [ -n "$out" ]; then
      c="${region:-CHROM}"
      { printf "$hdr"; printf "$c\t100\t.\tA\tG\t99\tPASS\t.\tGT\t0/1\n"; } \
        | case "$out" in *.gz) gzip -c ;; *) cat ;; esac > "$out"
    fi ;;
  index)     : > "$last.tbi" ;;
  --version) echo "bcftools 1.99-fake" ;;
esac
exit 0
"""


# Scores only contigs CADD's GRCh38 data actually knows -- Ensembl names, so MT and not M --
# and silently drops the rest, which is precisely what the real tool does: an unrecognised
# contig produces a fragment with a header and no rows, indistinguishable from an interval
# that genuinely had no variants.
#
# Echoing the input back instead would make this test VACUOUS. The old blunt mapping
# (`sed 's/^chr//'` in, `sed 's/^/chr/'` out) round-trips chrM -> M -> chrM perfectly well;
# what it got wrong was the name it handed to CADD in the middle. Only a fake that refuses
# to score `M` can see that.
CADD_ENSEMBL = r"""
prev=""; out=""
for a in "$@"; do [ "$prev" = "-o" ] && out="$a"; prev="$a"; done
for a in "$@"; do in="$a"; done
{
  printf '#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n'
  grep -v '^#' "$in" \
    | awk 'BEGIN{FS=OFS="\t"}
           $1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/ {print $1, $2, $4, $5, "0.8", "9.566"}'
} | gzip -c > "$out"
exit 0
"""


def merged_cadd_contigs(tmp, tag):
    """The set of contig names in the published 01_cadd table, or None if it is absent."""
    path = tmp / f"results_{tag}" / "01_cadd" / "caddv1.6.out.tsv.gz"
    if not path.is_file():
        return None
    with gzip.open(path, "rt") as fh:
        return {l.split("\t")[0] for l in fh if l.strip() and not l.startswith("#")}


def run_case(tmp, chroms, tag):
    """Runs the pipeline over `chroms` with CADD scoring on and everything else skipped."""
    binf = build_env(tmp)
    write_exec(binf / "bcftools", bcftools_fake(chroms))
    write_exec(tmp / "cadd_dir" / "CADD.sh", CADD_ENSEMBL)
    # SpliceAI/SQUIRLS skipped so the run COMPLETES (they need java and pysam, which no PATH
    # shim can supply) -- CADD itself is deliberately not skipped, since it is what is on test.
    trace, _named, out = run_pipeline(tmp, binf, ["--skip_spliceai_squirls"], tag)
    return trace, out


def main():
    if not shutil.which("nextflow"):
        print("SKIP: nextflow not on PATH")
        return 0

    # --- GATK-named callset: chr1/chrM must survive the round trip -------------------------
    print("chr-prefixed callset -- chrM must come back as chrM, not M or chrMT")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        trace, out = run_case(tmp, ["chr1", "chrM"], "chr")

        check("run_cadd and merge_cadd both ran",
              {"run_cadd", "merge_cadd"} <= trace,
              f"trace={sorted(trace)} -- nothing below proves anything if CADD did not run")
        check("the run completed",
              "[FAILED]" not in out and "ERROR" not in out, f"...{out.strip()[-400:]}")

        contigs = merged_cadd_contigs(tmp, "chr")
        check("a merged CADD table is published", contigs is not None, "no 01_cadd output")
        if contigs:
            check("the autosome round-trips", "chr1" in contigs, str(sorted(contigs)))
            check("the mitochondrion comes back as chrM",
                  "chrM" in contigs,
                  f"{sorted(contigs)} -- chrM was never scored, because run_cadd sent CADD "
                  "a contig name it does not know")
            check("no chrMT (a blanket prefix re-add)",
                  "chrMT" not in contigs, str(sorted(contigs)))
            check("no bare M or MT (the prefix was never restored)",
                  not ({"M", "MT"} & contigs), str(sorted(contigs)))

    # --- Ensembl-named callset: nothing may be rewritten ------------------------------------
    # chr_prefix is "" here, which is this repo's divergence from upstream. Upstream would
    # prefix every contig with chr and produce a table that joins nothing at all.
    print("Ensembl-named callset -- contigs must be returned untouched")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        trace, out = run_case(tmp, ["1", "MT"], "nochr")

        check("run_cadd and merge_cadd both ran",
              {"run_cadd", "merge_cadd"} <= trace, f"trace={sorted(trace)}")
        contigs = merged_cadd_contigs(tmp, "nochr")
        check("a merged CADD table is published", contigs is not None, "no 01_cadd output")
        if contigs:
            check("the autosome is not prefixed", "1" in contigs, str(sorted(contigs)))
            check("MT is left as MT", "MT" in contigs, str(sorted(contigs)))
            check("nothing gained a chr prefix",
                  not any(c.startswith("chr") for c in contigs), str(sorted(contigs)))

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("CADD contig names round-trip in both naming conventions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
