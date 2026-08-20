# Promote axis-A merge provenance to typed INFO keys, then strip the raw ones.
#
# Usage:  awk -f promote_caller_support.awk axisA.vcf axisA.vcf
#
# The file is named TWICE on purpose: pass 1 collects the contigs the records use so the
# header can declare the ones SVDB dropped (see the NR == FNR block). Naming it once prints
# nothing at all.
#
# This is mandatory, not tidiness. Verified by experiment during design: axis B carries
# axis A's INFO keys through VERBATIM and then appends its own, so a twice-merged record
# ends up with two svdb_origin, two FOUNDBY, two SUPP_VEC, two set and two VARID -- and
# bcftools returns the FIRST, which is the stale axis-A value. The same field then has
# different bit widths on different records of one file and every consumer silently reads
# the wrong one. tests/test_svdb_assumptions.py pins that upstream behaviour; this strips
# the collision before it can happen.
#
# What moves where:
#   svdb_origin=manta|delly  ->  CALLER_SUPP=manta,delly   which callers supported it
#   FOUNDBY=2                ->  NCALLER=2                 how many
#
# Records that never went through axis A already carry CALLER_SUPP/NCALLER from
# stamp_records.awk, so the field means one thing on every record of the cohort VCF. Those
# are left alone here -- no svdb_origin means nothing to promote.
#
# Axis B's OWN svdb_origin/SUPP_VEC are deliberately not touched by this: they name sample
# sets rather than callers, which is worth keeping and worth not confusing with caller
# support. Note SUPP_VEC's bits are ordered alphabetically by tag, not by --vcf order
# (verified in the design spike), so never read it positionally.

# MODE="cohort" switches this from axis-A behaviour to axis-B behaviour. They are not the
# same operation:
#
#   axis A  svdb_origin names CALLERS, so it promotes straight to CALLER_SUPP.
#   axis B  svdb_origin names SAMPLE SETS. Promoting it would overwrite caller support with
#           sample-set labels -- exactly the conflation warned about above, and what the first
#           run of this actually produced (CALLER_SUPP=COHORT,SAMP1,SAMP3). At axis B the
#           caller union is instead recovered from the per-input <tag>_INFO blobs, each of
#           which carries that input's own CALLER_SUPP, and axis B's svdb_origin/SUPP_VEC
#           are left alone as sample-set provenance worth keeping.
#
# SVDB rewrites the comma in a nested value as a colon, so a two-caller axis-A support
# arrives here as "CALLER_SUPP:tiddit:cnvpytor". Splitting on both is safe because
# check_caller_token() in samplesheet.nf rejects a caller containing either character.

BEGIN { FS = OFS = "\t" }

# Pass 1 -- the VCF is passed TWICE, so this reads it once to learn which contigs the records
# actually use and which the header already declares.
#
# SVDB copies one input's header into the merged VCF, so a contig that only the other input
# called on is never declared. `bcftools sort` does not survive that: it flushes a temp chunk
# written against the header it was handed, htslib then auto-adds the contig mid-stream, the
# two disagree, and the read dies with "Error encountered while parsing the input at chr19:1"
# -- a real failure on real data, where cnvnator called chr19 and the header SVDB kept did
# not name it. The give-away is those records sorting after every declared contig.
NR == FNR {
    if ($0 ~ /^##contig=<ID=/) {
        id = $0
        sub(/^##contig=<ID=/, "", id)
        sub(/[,>].*$/, "", id)
        declared_contig[id] = 1
    } else if ($0 !~ /^#/ && !($1 in seen_contig)) {
        seen_contig[$1] = 1
        contig_order[++n_contig] = $1
    }
    next
}

# Declare the promoted keys once, immediately before the column header.
/^#CHROM/ {
    # Backfill the contigs pass 1 found records on but the header never declared. Appended
    # rather than interleaved, so an undeclared contig sorts after the declared ones -- valid,
    # indexable, and not worth reordering a header to avoid.
    for (ci = 1; ci <= n_contig; ci++) {
        if (!(contig_order[ci] in declared_contig)) {
            print "##contig=<ID=" contig_order[ci] ">"
        }
    }
    if (!have_supp) {
        print "##INFO=<ID=CALLER_SUPP,Number=.,Type=String,Description=\"Callers supporting this record\">"
    }
    if (!have_ncaller) {
        print "##INFO=<ID=NCALLER,Number=1,Type=Integer,Description=\"Number of supporting callers\">"
    }
    print
    next
}

/^##INFO=<ID=CALLER_SUPP[,=]/ { have_supp = 1; print; next }
/^##INFO=<ID=NCALLER[,=]/     { have_ncaller = 1; print; next }
/^#/ { print; next }

{
    info = $8
    origin = ""
    foundby = ""
    kept = ""
    supp_value = ""
    cohort_supp = ""
    cohort_n = 0

    # SVDB's per-input blobs (<tag>_INFO, <tag>_SAMPLE, <tag>_FILTERS) copy a source
    # record's INFO, sample and FILTER columns into an INFO *value*, and it does not escape
    # every semicolon it copies. INFO is semicolon-delimited, so one such value arrives here
    # as two fields and the second is not KEY=VALUE at all -- bcftools rejects the record
    # ("not defined in the header") the moment it reads this back, which is what killed the
    # first real cohort run. Reattach the orphan to the field it broke off from, using
    # SVDB's own '|' separator.
    #
    # A single-token orphan -- what a FILTER of "MinQUAL;LowDepth" leaves behind -- is
    # indistinguishable from a legitimate INFO flag like IMPRECISE on its own, so it counts
    # as an orphan only when it follows a blob. Flags anywhere else are left alone.
    #
    # ponytail: a real flag written immediately after a blob is absorbed into it. SVDB
    # appends the blobs after the source record's own INFO, so that ordering does not arise
    # in its output; if some other producer ever emits one, the flag is lost, not corrupted.
    n = split(info, raw, ";")
    n_parts = 0
    for (i = 1; i <= n; i++) {
        is_field = (raw[i] ~ /^[A-Za-z_][A-Za-z0-9_.]*=/)
        is_flag = (raw[i] ~ /^[A-Za-z_][A-Za-z0-9_.]*$/)
        after_blob = (n_parts > 0 && parts[n_parts] ~ \
            /^[A-Za-z_][A-Za-z0-9_.]*_(CHROM|POS|QUAL|FILTERS|INFO|SAMPLE|FORMAT)=/)
        if (n_parts > 0 && !is_field && (!is_flag || after_blob)) {
            parts[n_parts] = parts[n_parts] "|" raw[i]
        } else {
            parts[++n_parts] = raw[i]
        }
    }
    n = n_parts

    # First pass: find the merge provenance, so the second pass knows whether an existing
    # CALLER_SUPP is stale (this record was merged) or authoritative (it was not).
    for (i = 1; i <= n; i++) {
        if (parts[i] ~ /^svdb_origin=/) origin = substr(parts[i], 13)
        else if (parts[i] ~ /^FOUNDBY=/) foundby = substr(parts[i], 9)
    }

    # At axis B the callers come from the per-input blobs, not from svdb_origin.
    if (MODE == "cohort") {
        delete seen
        union = ""
        ncall = 0
        for (i = 1; i <= n; i++) {
            if (parts[i] !~ /_INFO=/) continue
            m = split(parts[i], blob, "\\|")
            for (j = 1; j <= m; j++) {
                if (blob[j] !~ /^CALLER_SUPP:/) continue
                c = split(substr(blob[j], 13), callers, /[:,]/)
                for (k = 1; k <= c; k++) {
                    name = callers[k]
                    if (name == "" || (name in seen)) continue
                    seen[name] = 1
                    ncall += 1
                    union = (union == "") ? name : union "," name
                }
            }
        }
        origin = ""          # do not promote sample-set labels
        if (union != "") { cohort_supp = union; cohort_n = ncall }
        else { cohort_supp = ""; cohort_n = 0 }
    }

    for (i = 1; i <= n; i++) {
        kv = parts[i]
        if (kv == "" || kv == ".") continue

        # Raw merge keys: at axis A they are dropped so axis B's versions are the only ones
        # of that name. At axis B they are kept -- they name sample sets, which is
        # real provenance, and there is no later merge for them to collide with.
        if (MODE != "cohort") {
            if (kv ~ /^svdb_origin=/) continue
            if (kv ~ /^FOUNDBY=/)     continue
            if (kv ~ /^SUPP_VEC=/)    continue
            if (kv ~ /^set=/)         continue
            if (kv ~ /^VARID=/)       continue
        }

        # Superseded by what this merge actually found.
        if ((origin != "" || cohort_supp != "") && kv ~ /^CALLER_SUPP=/) continue
        if ((origin != "" || cohort_supp != "") && kv ~ /^NCALLER=/)     continue

        kept = (kept == "") ? kv : kept ";" kv
    }

    # No svdb_origin means this record did not come through a merge; its CALLER_SUPP and
    # NCALLER from standardization are already right and were kept above.
    if (origin != "") {
        gsub(/\|/, ",", origin)
        if (foundby == "") foundby = split(origin, tmp, ",")
        kept = (kept == "") ? "" : kept ";"
        kept = kept "CALLER_SUPP=" origin ";NCALLER=" foundby
        supp_value = origin
    } else if (cohort_supp != "") {
        kept = (kept == "") ? "" : kept ";"
        kept = kept "CALLER_SUPP=" cohort_supp ";NCALLER=" cohort_n
        supp_value = cohort_supp
    }

    # ALGORITHMS mirrors CALLER_SUPP after every merge. Stamped pre-merge it names one
    # caller, and SVDB keeps only the priority record's value -- so left alone it reports
    # one caller where three agreed, and the "ALGORITHMS populated on every record" check
    # passes while the field is wrong. Re-derived here it names everything that contributed.
    if (supp_value != "") {
        gsub(/(^|;)ALGORITHMS=[^;]*/, "", kept)
        sub(/^;/, "", kept)
        kept = (kept == "") ? "ALGORITHMS=" supp_value : kept ";ALGORITHMS=" supp_value
    }

    $8 = (kept == "") ? "." : kept
    print
}
