# List SVDB's per-input blob INFO keys from a VCF header as a bcftools `-x` list.
#
# Usage:  bcftools view -h cohort.vcf.gz | awk -f list_blob_info_keys.awk
#
# SVDB merge copies each input's CHROM/POS/QUAL/FILTER/INFO/FORMAT/SAMPLE columns into INFO
# values named after that input's tag: <tag>_INFO, <tag>_SAMPLE, <tag>_FILTERS and friends.
# AnnotSV copies INFO into the TSV verbatim, so an N-input cohort drops N whole source
# records into every row -- the bulk of an unreadable report, and enough on its own to push a
# cell past Excel's 32767-character limit.
#
# The tag is the input's own label, so these keys cannot be named in a param the way
# --annotsv_drop_info names PRPOS/PREND. They are read off the header instead.
#
# The suffix set is the one promote_caller_support.awk uses to recognise a blob -- the two
# have to stay in step, because that script is what harvests CALLER_SUPP out of these blobs
# before anything downstream is allowed to lose them.
#
# Prints nothing when the header declares no blobs, which the caller reads as "nothing to
# strip" -- a single-input cohort never went through a merge and has none.

# INFO header lines are ##INFO=<ID=KEY,Number=...>, so field 2 under this FS is "ID=KEY".
# Reading only field 2 keeps a Description containing '<' or ',' harmless.
BEGIN { FS = "[<,>]" }

/^##INFO=<ID=/ {
    id = $2
    sub(/^ID=/, "", id)
    if (id ~ /_(CHROM|POS|QUAL|FILTERS|INFO|SAMPLE|FORMAT)$/) {
        printf "%sINFO/%s", (n++ ? "," : ""), id
    }
}

END { if (n) print "" }
