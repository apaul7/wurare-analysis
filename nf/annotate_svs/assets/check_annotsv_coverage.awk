# Report which cohort VCF records reached the AnnotSV TSV, and which did not.
#
# Usage:  awk -f check_annotsv_coverage.awk cohort.vcf annotsv.tsv > coverage.txt
#
# The original acceptance criterion asked for "AnnotSV TSV row count matches the cohort VCF
# record count". That check is wrong in both directions and would fail on a perfectly correct
# run:
#
#   * The pipeline runs -annotationMode both, which emits full AND split rows -- one per SV
#     *and* one per SV x gene -- so TSV rows outnumber cohort records by design. Counting
#     DISTINCT IDs rather than rows is what makes this report survive that.
#   * Even in full mode the counts are only approximately equal. AnnotSV drops records below
#     its own -SVminSize and records it cannot type, so a legitimate run loses some.
#
# What matters is not the count but *which* records went missing, so this lists them by ID
# rather than asserting a number. Records legitimately drop out; silently losing one does
# not. The report is published, not merely printed.
#
# Exit status is always 0: a dropped record is information for the reader, not a pipeline
# failure. The one thing treated as an error is a TSV with no recognisable ID column, since
# that means the comparison never actually happened.

BEGIN { FS = OFS = "\t"; id_col = 0 }

# --- first file: the cohort VCF ---------------------------------------------------
FILENAME == ARGV[1] {
    if ($0 ~ /^#/ || $0 ~ /^[[:space:]]*$/) next   # header, or a stray blank line
    vcf_ids[$3] = 1
    n_vcf += 1
    next
}

# --- second file: the AnnotSV TSV -------------------------------------------------
FNR == 1 {
    # AnnotSV carries the input VCF's ID in a column named "ID". Located by name because
    # the column order differs between AnnotSV releases.
    for (i = 1; i <= NF; i++) if ($i == "ID") id_col = i
    if (id_col == 0) {
        print "ERROR: no 'ID' column in the AnnotSV TSV -- cannot verify coverage" > "/dev/stderr"
        exit 1
    }
    next
}

{
    if ($id_col in vcf_ids) { seen[$id_col] = 1 }
    else { unexpected[$id_col] = 1; n_unexpected += 1 }
    n_tsv += 1
}

END {
    if (id_col == 0) exit 1

    n_seen = 0
    for (id in seen) n_seen += 1

    print "cohort VCF records", n_vcf + 0
    print "AnnotSV TSV rows", n_tsv + 0
    print "cohort records present in TSV", n_seen

    n_missing = 0
    for (id in vcf_ids) if (!(id in seen)) { n_missing += 1; missing[n_missing] = id }
    print "cohort records absent from TSV", n_missing

    if (n_missing > 0) {
        print ""
        print "# absent record IDs -- expected causes are AnnotSV's -SVminSize floor and"
        print "# records it cannot type. Anything else is worth investigating."
        for (i = 1; i <= n_missing; i++) print missing[i]
    }

    if (n_unexpected > 0) {
        print ""
        print "# TSV rows whose ID is not in the cohort VCF -- should be none;"
        print "# non-zero means the TSV and the VCF are not from the same run."
        for (id in unexpected) print id
    }
}
