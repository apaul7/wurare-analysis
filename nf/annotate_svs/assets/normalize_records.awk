# Normalize the SV fields every merger assumes are uniform, and that callers disagree on.
#
# Usage:  awk -f normalize_records.awk in.vcf
#
# This replaces `svtk standardize`, which cannot run: its bundled no_contigs_template.vcf
# has a #CHROM line with no sample column, modern pysam rejects it as invalid VCF, and it is
# loaded unconditionally with no override. Broken in every current biocontainer build.
#
# What that means here: this is no longer a fallback branch for exotic callers, it is the
# only standardization the pipeline has, and it runs over manta and delly too. The three
# repairs below are the known silent-corruption sources.
#
#   SVLEN sign      DEL negative, everything else positive. Callers disagree, and size
#                   similarity matching then fails, so the same deletion lands in two
#                   clusters and is reported twice.
#   END missing     An absent or POS-equal END makes the event 1 bp to every downstream
#                   tool -- it never matches anything and never merges. Recomputed from
#                   SVLEN where SVLEN is trustworthy.
#   SVTYPE missing  Derived from a symbolic ALT. Without it a record is untyped and every
#                   type-aware merger drops or mis-clusters it.
#
# BND is left alone throughout: a breakend has no span, so END and SVLEN mean something
# different there and "fixing" them invents data. Whether BNDs stay in scope at all is the
# open decision -- this does not answer it, it just does not corrupt them.
#
# NOT done here, deliberately: DUP-encoded-as-INS reconciliation. Collapsing the two is
# a claim about biology, not a format repair -- the caller may genuinely mean an insertion.
# Jasmine's --dup_to_ins exists for it if it is ever wanted; doing it silently here would
# merge events that should stay distinct.

BEGIN { FS = OFS = "\t" }

/^#/ { print; next }

{
    svtype = ""; end = ""; svlen = ""
    n = split($8, parts, ";")
    for (i = 1; i <= n; i++) {
        if      (parts[i] ~ /^SVTYPE=/) svtype = substr(parts[i], 8)
        else if (parts[i] ~ /^END=/)    end    = substr(parts[i], 5)
        else if (parts[i] ~ /^SVLEN=/)  svlen  = substr(parts[i], 7)
    }

    # SVTYPE from a symbolic ALT when the caller omitted it. Digits and underscores are in the
    # class because real alleles carry them -- <CN0>/<CN2>/<CN3> from gCNV and GATK-SV,
    # <INS:ME:L1> from mobile-element callers. An alpha-only class typed none of them.
    if (svtype == "" && $5 ~ /^<[A-Z0-9:_]+>$/) {
        svtype = substr($5, 2, length($5) - 2)
        sub(/:.*/, "", svtype)          # <DUP:TANDEM> -> DUP
        # <CN#> are copy-number alleles, so the type is CNV. A literal "CN0" would be worse
        # than nothing -- SVDB and AnnotSV switch on SVTYPE and an unknown value takes their
        # default branch. NOT mapped to DEL/DUP: that asserts a direction the allele lacks.
        if (svtype ~ /^CN[0-9]+$/) svtype = "CNV"
    }

    # BND and INS are both exempt from the END/SVLEN derivations, for the same reason: neither
    # spans reference. END == POS is spec-correct for an insertion, not a collapsed 1 bp
    # event, so deriving END = POS + SVLEN invents an interval the variant never occupied --
    # which then clusters against real DELs and DUPs and drags in genes it does not touch.
    if (svtype != "" && svtype != "BND" && svtype != "INS") {
        abslen = (svlen + 0 < 0) ? -(svlen + 0) : (svlen + 0)

        # END from SVLEN when END is absent or collapsed to a 1 bp event.
        if ((end == "" || end + 0 <= $2 + 0) && abslen > 0) {
            end = $2 + abslen
        }
        # SVLEN from END when the caller gave only a span.
        if ((svlen == "" || abslen == 0) && end != "" && end + 0 > $2 + 0) {
            abslen = end - $2
        }
        # Sign by type: DEL negative, everything else positive.
        if (abslen > 0) {
            svlen = (svtype == "DEL") ? -abslen : abslen
        }
    }
    # An insertion still gets its sign normalized -- the length of inserted sequence is always
    # positive -- but keeps whatever END the caller wrote.
    else if (svtype == "INS" && svlen != "") {
        abslen = (svlen + 0 < 0) ? -(svlen + 0) : (svlen + 0)
        if (abslen > 0) svlen = abslen
    }

    # END before POS survives the repairs above only when there was no SVLEN to rebuild it
    # from. It is not ignorable: htslib rejects the record when the VCF is read back, in a
    # message naming neither the caller nor the record. Say which record, here, while the
    # input file still identifies the caller that produced it.
    if (end != "" && end + 0 < $2 + 0) {
        n_bad_end += 1
        if (n_bad_end <= 5) {
            printf "normalize_records: WARNING %s:%s (%s) has END=%s before POS and no SVLEN " \
                   "to repair it from\n", $1, $2, $3, end > "/dev/stderr"
        }
    }

    # Rebuild INFO with the normalized values, preserving every other key and its order.
    out = ""
    seen_type = 0; seen_end = 0; seen_len = 0
    for (i = 1; i <= n; i++) {
        kv = parts[i]
        if (kv == "" || kv == ".") continue
        if      (kv ~ /^SVTYPE=/) { if (svtype == "") continue; kv = "SVTYPE=" svtype; seen_type = 1 }
        else if (kv ~ /^END=/)    { if (end == "")    continue; kv = "END=" end;       seen_end  = 1 }
        else if (kv ~ /^SVLEN=/)  { if (svlen == "")  continue; kv = "SVLEN=" svlen;   seen_len  = 1 }
        out = (out == "") ? kv : out ";" kv
    }
    if (!seen_type && svtype != "") out = (out == "") ? "SVTYPE=" svtype : out ";SVTYPE=" svtype
    if (!seen_end  && end   != "")  out = (out == "") ? "END=" end       : out ";END=" end
    if (!seen_len  && svlen != "")  out = (out == "") ? "SVLEN=" svlen   : out ";SVLEN=" svlen

    $8 = (out == "") ? "." : out
    print
}

END {
    if (n_bad_end > 5) {
        printf "normalize_records: WARNING %d record(s) total have END before POS; the first " \
               "5 are named above\n", n_bad_end > "/dev/stderr"
    }
}
