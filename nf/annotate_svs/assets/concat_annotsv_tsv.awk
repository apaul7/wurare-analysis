# Concatenate per-shard AnnotSV TSVs into one, keeping exactly one header line.
#
# Usage:  awk -f concat_annotsv_tsv.awk shard_00001.annotsv.tsv shard_00002.annotsv.tsv ...
#
# The shards are the per-shard AnnotSV outputs produced after shard_vcf.awk split the cohort
# VCF to get under Tcl's 2 GiB value limit. See that file for why the split exists at all.
#
# THE HEADER MUST APPEAR ONCE, AND FIRST. This is not tidiness. Both consumers of the
# assembled TSV locate their columns BY NAME from line 1 -- check_annotsv_coverage.awk reads
# it at FNR==1 to find the ID column, filter_annotsv_tsv.awk at NR==1 to build its whole
# column map -- because AnnotSV's column set shifts between releases and with the annotation
# bundle. A header repeated mid-file is therefore not a cosmetic wart: it becomes a data row
# whose every cell is a column name, the tier filters evaluate it as a variant, and the
# coverage report lists "ID" among the records missing from the cohort VCF.
#
# EVERY SHARD'S HEADER MUST MATCH THE FIRST. Shards annotated with different AnnotSV versions
# or bundles would have different column sets, and concatenating those bodies under one header
# would shift columns silently from the first mismatched shard onward -- every downstream
# value read from the wrong column, with nothing to show for it. Cheap to check, so checked.
#
# Empty and header-only shards contribute nothing and are not an error: a shard whose records
# AnnotSV could not type legitimately yields no rows.

FNR == 1 {
    if (!seen) {
        header = $0
        print
        seen = 1
    }
    else if ($0 != header) {
        printf("ERROR: %s has a different header from the first shard.\n", FILENAME) > "/dev/stderr"
        printf("  first shard: %s\n", header) > "/dev/stderr"
        printf("  this shard:  %s\n", $0) > "/dev/stderr"
        printf("  concatenating these would shift every column from here on, silently.\n") > "/dev/stderr"
        bad = 1
        exit 1
    }
    next
}

{ print }

END {
    if (bad) exit 1
    # No header at all means every shard was empty, which downstream reads as "no ID column"
    # several steps later. Fail here, where the cause is still visible.
    if (!seen) {
        print "ERROR: no shard produced a header -- nothing to concatenate" > "/dev/stderr"
        exit 1
    }
}
