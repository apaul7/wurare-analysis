# The bcftools `annotate -x` keep-list for the Talos handoff VCF.
#
# Usage:  awk -v GNOMAD_POP=gnomad_v4.1 -f talos_keep_info.awk cohort.schema.vcf
#
# Prints one line: ^INFO/TAG,...,^FORMAT/TAG,... -- the ^-prefixed KEEP form of
# `bcftools annotate -x`, which removes every tag of that class not listed, values and
# header lines both. Restricted to tags the header actually declares, because bcftools
# errors on a tag the file lacks and e.g. CHR2 is only present when a BND survived this
# far. FILTER is deliberately untouched: it is already PASS after talos_schema.awk.
#
# FORMAT is stripped to GT plus the duphold depths -- everything Talos reads -- and that
# strip is load-bearing, not cosmetic. svdb --merge keeps ONE input's ##FORMAT
# declarations while records keep their own caller's sample text, and callers disagree
# on the same ID: manta writes SR/PR as Number=. ref,alt pairs ("98,0"), lumpy declares
# SR as Number=1 Integer. bcftools shrugs at the mismatch; hail dies loading the handoff
# ("parseFormatInt: invalid character ','"). Dropping the caller FORMAT baggage removes
# the lying declarations and the values they lie about in one move.
#
# Why strip at all: the cohort VCF accumulates every input caller's INFO baggage --
# PRPOS/PREND probability arrays, caller-private keys -- and all of it rides into
# 06_talos otherwise. Extra fields are not errors for Talos (hail reads the header and
# ignores what nothing accesses); this is for the human reading the handoff file and for
# its size, not for correctness.
#
# Keep set = the fields check_talos_fields.awk audits (Talos's read set, traced from its
# rearrange_annotations()), plus SOFT_FILTERS (this pipeline's own carry), plus
# gnomad_sv_AF (the documented source of the {GNOMAD_POP}_sv_AF copy -- talos_schema.awk
# explains the dual name), plus every PREDICTED_* tag gatk SVAnnotate wrote, wholesale,
# because Talos's consequence logic reads that family rather than PREDICTED_LOF alone.
# FORMAT keep set = GT plus the duphold depths merge_depth.awk declared, which the
# Talos report's evidence panel reads.
#
# Runs on talos_schema.awk output, which always declares SVTYPE and friends -- so an
# empty keep-list means the input was not that output, and exiting nonzero with no
# stdout is what makes the caller's `KEEP=$(awk ...)` fail under `set -e` instead of
# handing bcftools a bare "^" that would strip everything.

BEGIN {
    if (GNOMAD_POP == "") GNOMAD_POP = "gnomad_v4.1"
    split("SVTYPE SVLEN END CHR2 END2 ALGORITHMS STATUS AC AF AN " \
          "AF_MALE AF_FEMALE N_HET N_HOMALT SOFT_FILTERS gnomad_sv_AF", w, " ")
    for (i in w) wanted[w[i]] = 1
    wanted[GNOMAD_POP "_sv_AF"] = 1
    wanted[GNOMAD_POP "_sv_SVID"] = 1
    n = 0
    split("GT DHFFC DHBFC DHFC DHBZ", fw, " ")
    for (i in fw) fwanted[fw[i]] = 1
    fn = 0
}

/^##INFO=<ID=/ {
    id = $0
    sub(/^##INFO=<ID=/, "", id)
    sub(/[,>].*$/, "", id)
    if (id in wanted || id ~ /^PREDICTED_/) keep[n++] = id
    next
}

/^##FORMAT=<ID=/ {
    id = $0
    sub(/^##FORMAT=<ID=/, "", id)
    sub(/[,>].*$/, "", id)
    if (id in fwanted) fkeep[fn++] = id
    next
}

# Header is over; records carry no declarations.
/^#CHROM/ { exit }

END {
    if (n == 0) exit 1
    out = ""
    for (i = 0; i < n; i++) out = out (i ? "," : "") "INFO/" keep[i]
    # A file with no FORMAT declarations gets no FORMAT clause -- a bare "^FORMAT"
    # would be as destructive as the bare "^" the INFO guard exists for.
    for (i = 0; i < fn; i++) out = out "," (i ? "" : "^") "FORMAT/" fkeep[i]
    print "^" out
}
