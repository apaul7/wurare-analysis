# Report whether a cohort VCF carries every field Talos requires, before handing it over.
#
# Usage:  awk -v GNOMAD_POP=gnomad_v4.1 -f check_talos_fields.awk cohort.vcf
#
# Talos's failure mode is silence, which is why this exists. Two examples from
# tracing Talos's ingest code: a variant with no PREDICTED_LOF entry is dropped entirely
# rather than reported, and an absent ALGORITHMS is DEFAULTED to ['gCNV'] rather than
# rejected -- so a callset can pass through Talos and come out mislabelled or half empty
# with nothing anywhere saying so.
#
# Reports rather than fails: which fields are present, which are missing, and how many
# records carry PREDICTED_LOF (the hard gate). Exit 0 always -- the answer is for a human
# deciding whether the handoff is ready, and a pipeline that refused to publish the report
# would be withholding exactly the diagnostic that is wanted.
#
# Two questions, not one, because Talos fails two different ways:
#
#   HEADER   hail builds mt.info from the VCF header, and Talos reads most fields as a direct
#            struct access. An UNDECLARED field is an error before any filtering happens --
#            Talos does not start. Only ALGORITHMS, STATUS, CHR2 and END2 are tolerated
#            absent, per its rearrange_annotations().
#   RECORDS  a declared field with no value on a record reads as missing, which is the silent
#            path: PREDICTED_LOF absent drops that variant, ALGORITHMS absent defaults the
#            whole callset to 'gCNV'.
#
# A field can be declared and still useless, and an undeclared one is fatal regardless of what
# the records hold. Hence two tables answering two different questions.

BEGIN {
    FS = "\t"
    if (GNOMAD_POP == "") GNOMAD_POP = "gnomad_v4.1"
    split("SVTYPE,SVLEN,END,ALGORITHMS,STATUS,PREDICTED_LOF,AC,AF,AN,N_HET,N_HOMALT", req, ",")
    req[12] = GNOMAD_POP "_sv_AF"
    req[13] = GNOMAD_POP "_sv_SVID"
    nreq = 13

    # The fields Talos reads without first checking they exist. Undeclared, it raises.
    split("SVTYPE,SVLEN,END,PREDICTED_LOF,AC,AF,AN,N_HET,N_HOMALT", fatal, ",")
    fatal[10] = GNOMAD_POP "_sv_AF"
    fatal[11] = GNOMAD_POP "_sv_SVID"
    nfatal = 11
}

/^##INFO=<ID=/ {
    id = $0
    sub(/^##INFO=<ID=/, "", id)
    sub(/[,>].*$/, "", id)
    declared[id] = 1
    next
}

/^#/ { next }

{
    n_records += 1
    n = split($8, parts, ";")
    delete here
    for (i = 1; i <= n; i++) {
        eq = index(parts[i], "=")
        key = (eq > 0) ? substr(parts[i], 1, eq - 1) : parts[i]
        here[key] = 1
    }
    for (r = 1; r <= nreq; r++) if (req[r] in here) seen[req[r]] += 1
    if ("PREDICTED_LOF" in here) n_lof += 1
}

END {
    print "records\t" n_records + 0
    print ""
    print "# Header declaration. Talos reads these off mt.info directly, and hail builds that"
    print "# struct from the header -- an UNDECLARED field raises before filtering starts."
    print "field\theader"
    n_undeclared = 0
    for (r = 1; r <= nfatal; r++) {
        ok = (fatal[r] in declared)
        if (!ok) n_undeclared += 1
        print fatal[r] "\t" (ok ? "declared" : "UNDECLARED")
    }

    # Talos branches on AF_MALE and then reads AF_FEMALE off the same struct, so the pair has
    # to arrive together. Half a pair is the crash case -- and it is the likely one, because
    # bcftools +fill-tags declares only the sex groups the PED actually contains.
    sexed = (("AF_MALE" in declared) && ("AF_FEMALE" in declared)) ? "AF_MALE/AF_FEMALE" \
          : ((("MALE_AF" in declared) && ("FEMALE_AF" in declared)) ? "MALE_AF/FEMALE_AF" : "")
    if (sexed == "") n_undeclared += 1
    print "sex-stratified AF pair\t" (sexed == "" ? "UNDECLARED" : "declared as " sexed)

    if (n_undeclared > 0) {
        print ""
        print "# WARNING: " n_undeclared " field(s) undeclared -- Talos raises on load rather"
        print "# than filtering badly. Fix the header before handing this VCF over."
    }

    print ""
    print "field\trecords_with_it\tstatus"
    for (r = 1; r <= nreq; r++) {
        c = seen[req[r]] + 0
        status = (c == 0) ? "MISSING" : ((c < n_records) ? "partial" : "ok")
        print req[r] "\t" c "\t" status
    }
    print ""
    print "# PREDICTED_LOF is Talos's hard gate: a variant without at least one entry is"
    print "# dropped entirely, not reported. Only gatk SVAnnotate produces it."
    print "records_retainable_by_talos\t" n_lof + 0
    if (n_records > 0 && n_lof == 0) {
        print ""
        print "# WARNING: no record carries PREDICTED_LOF -- Talos would retain nothing."
    }
}
