# Rewrite sequence-resolved ALTs to symbolic ones, for gatk SVAnnotate only.
#
# Usage:  awk -f symbolic_alt.awk cohort.vcf
#
# SVAnnotate refuses any record whose ALT is literal sequence:
#
#   java.lang.IllegalArgumentException: Unexpected ALT allele: TTTTTTCTTTCTTT...
#   Expected breakpoint or symbolic ALT allele representing a structural variant record.
#
# and it refuses by dying, so a single such record takes the whole Talos tail with it. Manta
# emits these routinely -- an insertion whose sequence it assembled, or a deletion written out
# in full rather than as <DEL> -- so this is ordinary caller output, not malformed input.
#
# Talos-side only, deliberately. A literal ALT carries more than <DEL> does, and AnnotSV,
# duphold and every published VCF read it happily; GATK is the one consumer that cannot. So
# the rewrite happens on the branch feeding SVAnnotate and nowhere upstream, exactly as the
# FILTER-to-SOFT_FILTERS move does. 04_annotate and 05_filter keep the sequence.
#
# What that costs, stated rather than buried: the inserted or deleted bases are gone from
# 06_talos. Talos reads SVTYPE, SVLEN, END and PREDICTED_LOF and never the ALT sequence, so
# nothing downstream of here wants them -- but if a consumer ever does, they are in the files
# above.
#
# Untouched: symbolic ALTs (<DEL>, <DUP:TANDEM>), breakends (N[chr2:1000[), the missing ALTs
# "." and "*", and multi-allelic records. Two cases cannot be converted at all, and both pass
# through unchanged and counted to stderr rather than guessed at:
#
#   no SVTYPE     nothing names the symbol. Inferring DEL from REF/ALT lengths would be
#                 inventing a call.
#   multi-allelic a record carries ONE SVTYPE and there is no honest way to spread it over
#                 several alleles. Writing a single <SVTYPE> over the column would collapse
#                 distinct alleles into one -- wrong, and wrong quietly.
#
# SVAnnotate still rejects both, which is the right outcome for a record nobody can convert:
# it fails loudly rather than being silently mangled first.

function info_get(info, key,    n, i, parts) {
    n = split(info, parts, ";")
    for (i = 1; i <= n; i++) if (parts[i] ~ ("^" key "=")) return substr(parts[i], length(key) + 2)
    return ""
}

# One allele, not a whole ALT column: a symbolic allele, a breakend, and the two missing
# forms are all things SVAnnotate accepts as-is.
function is_literal(allele) {
    return (allele !~ /^</ && allele !~ /[\[\]]/ && allele != "." && allele != "*")
}

function any_literal(alt,    n, i, alleles) {
    n = split(alt, alleles, ",")
    for (i = 1; i <= n; i++) if (is_literal(alleles[i])) return 1
    return 0
}

BEGIN { FS = OFS = "\t" }

/^#/ { print; next }

# Multi-allelic records are left alone, because a record carries ONE SVTYPE and there is no
# honest way to spread it over several alleles. Rewriting the column to a single <SVTYPE>
# would collapse two distinct alleles into one -- wrong, and wrong quietly, which is the
# failure class this pipeline exists to avoid. So the record passes through and, if any of
# its alleles is literal, is counted for the same warning as an untyped one: SVAnnotate will
# still reject it, loudly, which is the correct outcome for something nobody can convert.
$5 ~ /,/ {
    if (any_literal($5)) multiallelic += 1
    print
    next
}

{
    alt = $5

    if (is_literal(alt)) {
        svtype = info_get($8, "SVTYPE")
        if (svtype != "") {
            $5 = "<" svtype ">"
            # A symbolic ALT means the REF is the single anchoring base at POS, not the whole
            # deleted span. A long REF beside a symbolic ALT is a record stating two different
            # things about its own length.
            $4 = substr($4, 1, 1)
            converted += 1
        }
        else {
            untyped += 1
        }
    }

    print
}

END {
    if (converted > 0) {
        printf "symbolic_alt: rewrote %d sequence-resolved ALT(s) to symbolic\n", \
            converted > "/dev/stderr"
    }
    if (untyped > 0) {
        printf "symbolic_alt: WARNING %d record(s) have a literal ALT and no SVTYPE; " \
               "SVAnnotate will reject them\n", untyped > "/dev/stderr"
    }
    if (multiallelic > 0) {
        printf "symbolic_alt: WARNING %d multi-allelic record(s) carry a literal ALT; " \
               "one SVTYPE cannot describe several alleles, so they are left unchanged and " \
               "SVAnnotate will reject them\n", multiallelic > "/dev/stderr"
    }
}
