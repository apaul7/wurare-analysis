# Stamp a unique record ID and an authoritative INFO/ALGORITHMS onto a VCF.
#
# Usage:  awk -v TAG=<sample_set>.<caller> -v CALLER=<caller> -f stamp_records.awk in.vcf
#
# Lives in its own file rather than inline in the process script for one reason: an awk
# program embedded in a Nextflow triple-quoted string needs every $ and " escaped past
# Groovy interpolation and then past the shell, which is how the first version of this
# silently broke. Here it is plain awk, and tests/test_prepare_stamps.py runs it directly
# with no Nextflow and no container.
#
# ID: a running counter, not CHROM_POS -- two callers can and do emit two records at one
# position, and the ID has to survive as a join key for duphold's per-sample FORMAT
# recombination and for `bcftools merge -m id` on the truvari upgrade path.
#
# ALGORITHMS: the samplesheet's caller wins, always. An input that already carries an
# ALGORITHMS is restamped, not deferred to -- inputs do arrive pre-stamped (a re-run, or a
# VCF that passed through another pipeline), and honouring a stale value silently
# mislabels the whole callset. The design requires it be populated from the caller column and
# never defaulted; Talos defaults absent ALGORITHMS to gCNV, so wrong here is not rejected
# downstream, just wrong.

BEGIN { FS = OFS = "\t"; n = 0; have_header = 0; have_supp = 0; have_ncaller = 0 }

# Remember an existing declaration so the header is not emitted twice.
/^##INFO=<ID=ALGORITHMS[,=]/  { have_header = 1; print; next }
/^##INFO=<ID=CALLER_SUPP[,=]/ { have_supp = 1; print; next }
/^##INFO=<ID=NCALLER[,=]/     { have_ncaller = 1; print; next }

/^##/ { print; next }

/^#CHROM/ {
    if (!have_header) {
        print "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Source algorithms\">"
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

{
    n += 1
    $3 = TAG "_" n

    info = $8
    # Strip any existing ALGORITHMS, wherever it sits in the INFO string, then re-add.
    gsub(/(^|;)ALGORITHMS=[^;]*/, "", info)
    sub(/^;/, "", info)
    if (info == "." || info == "") info = "ALGORITHMS=" CALLER
    else info = info ";ALGORITHMS=" CALLER

    # CALLER_SUPP/NCALLER are stamped here as well as promoted after a merge, so that a
    # record which never goes through axis A -- a joint VCF entering axis B directly, or a
    # caller with no partner on the same sample set -- still carries them. Without this the
    # field is present on some cohort records and absent on others, and a filter written
    # against it means different things on different rows.
    gsub(/(^|;)CALLER_SUPP=[^;]*/, "", info)
    gsub(/(^|;)NCALLER=[^;]*/, "", info)
    sub(/^;/, "", info)
    info = info ";CALLER_SUPP=" CALLER ";NCALLER=1"
    sub(/^;/, "", info)

    $8 = info

    print
}
