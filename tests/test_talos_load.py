#!/usr/bin/env python3
"""Does Talos actually load what the tail emits?

Every other Talos check in this repo reasons about the schema from the outside --
check_talos_fields.awk reads a VCF and says what it thinks Talos will make of it.
This one asks Talos. It imports the real `rearrange_annotations` out of the pinned
image and runs it over a VCF that has been through talos_schema.awk, which is the
exact call that raised before the header backfill existed:

    KeyError: "StructExpression instance has no field 'gnomad_v4.1_sv_SVID'"

Talos publishes no image, so build theirs first (their README: docker build -t talos:<v> .):

    TALOS_IMAGE=talos:11.1.0 python3 tests/test_talos_load.py

Skips, loudly, when the image is absent -- a check nobody can build that fails anyway is
worse than one that says why it did nothing. amd64 only: hail ships no aarch64 wheel, so
this runs under emulation on an arm64 Mac and takes a couple of minutes.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ASSETS = Path(__file__).parent.parent / "nf" / "annotate_svs" / "assets"
SCHEMA = ASSETS / "talos_schema.awk"
IMAGE = os.environ.get("TALOS_IMAGE", "talos:11.1.0")

# A cohort record carrying what the pipeline puts on it by the time the tail runs: AC/AN/AF
# from fill_tags at Phase 2, ALGORITHMS from Phase 1, PREDICTED_LOF from gatk SVAnnotate,
# AF_MALE/AF_FEMALE from fill-tags -S, gnomad_sv_AF from svdb query. Deliberately NOT
# carrying gnomad_v4.1_sv_SVID or the Talos-side names -- that is what the awk is for.
COHORT_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">
##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##INFO=<ID=AF_MALE,Number=A,Type=Float,Description="AF in MALE">
##INFO=<ID=AF_FEMALE,Number=A,Type=Float,Description="AF in FEMALE">
##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Algorithms">
##INFO=<ID=PREDICTED_LOF,Number=.,Type=String,Description="LoF genes">
##INFO=<ID=gnomad_sv_AF,Number=1,Type=Float,Description="svdb query AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP1\tSAMP2\tSAMP3
chr1\t1000\tr1\tN\t<DEL>\t99\tPASS\tSVTYPE=DEL;SVLEN=-100;END=1100;AC=1;AN=6;AF=0.17;\
AF_MALE=0.25;AF_FEMALE=0;ALGORITHMS=manta;PREDICTED_LOF=GENE1;gnomad_sv_AF=0.001\tGT\t0/1\t0/0\t0/0
"""

# Talos's own entry point needs a config, a PED and PanelApp data before it will run. The
# schema contract does not: rearrange_annotations is where every INFO field is read, and
# _force_count is what makes hail evaluate that struct rather than only plan it.
PROBE = """import sys
import hail as hl
from talos.run_hail_filtering_sv import rearrange_annotations

hl.init(master='local[1]', quiet=True)
mt = hl.import_vcf(sys.argv[1], force=True, reference_genome='GRCh38', skip_invalid_loci=True)
mt = rearrange_annotations(mt, hl.literal({'GENE1': 'ENSG1'}))
mt.rows()._force_count()
print('LOADED OK')
"""

failures = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def talos_loads(workdir, vcf_name):
    """Run Talos's annotation rearrangement over a VCF. Returns (ok, last error line)."""
    res = subprocess.run(
        ["docker", "run", "--rm", "--platform", "linux/amd64",
         "-v", f"{workdir}:/s", "-w", "/s", IMAGE, "python", "probe.py", vcf_name],
        capture_output=True, text=True)
    out = res.stdout + res.stderr
    if "LOADED OK" in out:
        return True, ""
    err = [l for l in out.splitlines() if "Error" in l or "error:" in l]
    return False, (err[-1][:200] if err else out.strip().splitlines()[-1][:200])


def main():
    if not SCHEMA.is_file():
        sys.exit(f"missing {SCHEMA}")
    if not shutil.which("docker"):
        print(f"SKIP: docker not on PATH -- cannot run {IMAGE}")
        return 0
    if subprocess.run(["docker", "image", "inspect", IMAGE],
                      capture_output=True).returncode != 0:
        print(f"SKIP: image {IMAGE} not present. Talos publishes none -- clone "
              "populationgenomics/talos and `docker build -t talos:<version> .`, or point "
              "TALOS_IMAGE at one you already have.")
        return 0

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "probe.py").write_text(PROBE)
        (td / "cohort.vcf").write_text(COHORT_VCF)
        awked = subprocess.run(
            ["awk", "-v", "GNOMAD_POP=gnomad_v4.1", "-f", str(SCHEMA), str(td / "cohort.vcf")],
            capture_output=True, text=True, check=True)
        (td / "cohort.talos.vcf").write_text(awked.stdout)

        print(f"the tail's output loads in {IMAGE}")
        ok, err = talos_loads(td, "cohort.talos.vcf")
        check("Talos reads every field it asks for", ok, err)

        # The regression this exists for. Talos reads most INFO fields as direct struct
        # accesses and hail builds that struct from the header, so the untransformed cohort
        # VCF -- valid and complete by every other measure -- raises before filtering starts.
        print("and the untransformed cohort VCF still does not")
        ok, err = talos_loads(td, "cohort.vcf")
        check("a VCF without the backfilled headers is rejected", not ok,
              "it loaded -- either Talos started tolerating absent fields, or the fixture "
              "already carries them, and either way the backfill needs rechecking")
        check("and it fails on the field the backfill adds",
              not ok and "gnomad_v4.1_sv_SVID" in err, f"got: {err}")

    print()
    if failures:
        print(f"{len(failures)} FAILED: {', '.join(failures)}")
        return 1
    print("Talos loads the tail's output")
    return 0


if __name__ == "__main__":
    sys.exit(main())
