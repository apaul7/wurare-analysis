# Split a VCF body into size-bounded shards, each a complete VCF with the header prepended.
#
# Usage:  bcftools view -H in.vcf.gz \
#           | awk -v MAXBYTES=1000000000 -v HDR=header.vcf -f shard_vcf.awk
# Emits:  shard_00001.vcf, shard_00002.vcf, ... in the working directory
#         and the shard COUNT on stdout, for the caller to assert against
#
# WHY THIS EXISTS. AnnotSV is written in Tcl, and Tcl 8 caps a single value at 2^31-1 bytes.
# A 300-sample cohort VCF is ~9.2 GB even after the SVDB blobs and --annotsv_drop_info keys
# are stripped, so AnnotSV dies in its very first step with "max size for a Tcl value
# (2147483647 bytes) exceeded". That is a ceiling, not a sizing problem: more memory cannot
# move it. Sharding the input is the way past it without dropping the genotype columns, which
# is what `bcftools view -G` would have cost.
#
# THE BUDGET IS BYTES, NOT RECORDS, because the limit is a byte limit and bytes-per-record is
# not a constant. It scales with sample count: a 300-sample cohort spends ~12.6 KB per record,
# nearly all of it genotype columns, where a handful of samples spends a few hundred bytes. A
# record-count knob tuned for one cohort is silently wrong for the next, and the failure lands
# hours into a run. This was learned the expensive way -- a 200000-record default, set from a
# GUESSED record count, produced a 2,520,953,496-byte shard: 17% over the ceiling.
#
# WHY NOT `split`. The bcftools container ships BusyBox split, which has no --filter, no -d
# and no --additional-suffix, so the header cannot be prepended per chunk that way. Checked
# rather than assumed.
#
# ORDER IS PRESERVED. Records are emitted in input order and every record lands in exactly
# one shard, so a coordinate-sorted input yields shards that are each a contiguous coordinate
# range. Contigs stay grouped and the concatenated report stays in coordinate order.

BEGIN {
    if (MAXBYTES == "" || MAXBYTES + 0 <= 0) {
        print "ERROR: shard_vcf.awk needs -v MAXBYTES=<max bytes per shard>" > "/dev/stderr"
        bad = 1
        exit 1
    }
    if (HDR == "") {
        print "ERROR: shard_vcf.awk needs -v HDR=<header file>" > "/dev/stderr"
        bad = 1
        exit 1
    }

    # The header is repeated in every shard, so it counts against every shard's budget. Read
    # once here: on a 300-sample cohort the #CHROM line alone is several KB.
    hdr_bytes = 0
    n_hdr = 0
    while ((r = (getline line < HDR)) > 0) {
        hdr_lines[++n_hdr] = line
        hdr_bytes += length(line) + 1
    }
    if (r < 0) {
        printf("ERROR: cannot read header file '%s'\n", HDR) > "/dev/stderr"
        bad = 1
        exit 1
    }
    if (n_hdr == 0) {
        printf("ERROR: header file '%s' is empty\n", HDR) > "/dev/stderr"
        bad = 1
        exit 1
    }
    close(HDR)

    if (hdr_bytes >= MAXBYTES + 0) {
        printf("ERROR: the VCF header alone is %d bytes, at or over the %d-byte budget\n",
               hdr_bytes, MAXBYTES) > "/dev/stderr"
        bad = 1
        exit 1
    }
}

{
    len = length($0) + 1

    # A single record that cannot fit in a shard of its own is unsplittable -- there is no
    # smaller unit than a VCF record here. Say so rather than emitting an over-budget shard
    # and letting AnnotSV fail on it later with a message that names neither.
    if (hdr_bytes + len > MAXBYTES + 0) {
        printf("ERROR: record %d is %d bytes, which does not fit a %d-byte shard alongside\n",
               NR, len, MAXBYTES) > "/dev/stderr"
        printf("       a %d-byte header. Raise the budget above %d.\n",
               hdr_bytes, hdr_bytes + len) > "/dev/stderr"
        bad = 1
        exit 1
    }

    if (out == "" || cur_bytes + len > MAXBYTES + 0) {
        if (out != "") close(out)
        idx += 1
        out = sprintf("shard_%05d.vcf", idx)
        for (i = 1; i <= n_hdr; i++) print hdr_lines[i] > out
        cur_bytes = hdr_bytes
    }

    print > out
    cur_bytes += len
}

END {
    if (bad) exit 1
    if (out != "") close(out)
    if (idx == 0) {
        print "ERROR: no records on stdin -- nothing to shard" > "/dev/stderr"
        exit 1
    }
    print idx
}
