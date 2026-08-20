# Recombine per-sample duphold FORMAT fields onto the cohort VCF, joined on record ID.
#
# Usage:  awk -f merge_depth.awk depth.tsv cohort.vcf
#
#   depth.tsv   ID <tab> SAMPLE <tab> DHFFC <tab> DHBFC <tab> DHFC <tab> DHBZ
#   cohort.vcf  the multi-sample VCF the per-sample extracts were cut from
#
# BOTH INPUTS MUST BE SORTED BY ID. See "the precondition" below -- it is checked here, not
# assumed, and a violation is fatal.
#
# The join is on ID and nothing else: duphold takes one
# alignment at a time, so a cohort VCF has to be split per sample and put back together,
# and a positional or POS/REF/ALT join is exactly how depth lands on the wrong record.
# Symbolic ALTs make that a real risk rather than a theoretical one -- two records can
# share CHROM/POS/REF/ALT and differ only in INFO/END, which bcftools merge will happily
# collapse. The unique ID stamped in Phase 1 (assets/stamp_records.awk) exists for this.
#
# A sample with no alignment gets "." for every field rather than being dropped: the
# failure mode to guard against is depth attached to the wrong sample, and a missing column silently
# shifting every later one is how that happens. Samples are matched by NAME, taken from the
# cohort VCF's own #CHROM line, so the order of the per-sample duphold runs cannot matter.
#
# --- why this is a merge join rather than a lookup table -----------------------------
# This used to load every (record, sample) pair into an awk array before reading the first
# VCF line. That is O(records x samples), and it does not survive a real cohort: at 300
# samples and ~3M records the depth table was 72 GB and ~900 million rows, needing an
# estimated 150 GB of RAM for the array alone. It did not fail -- it thrashed for 15 hours
# and produced nothing. Only the rows belonging to the record currently being written are
# held here, so peak memory is O(samples): a few hundred entries.
#
# --- the precondition, and why it is enforced ---------------------------------------
# A merge join needs both sides in the same order. recombine_depth sorts them, under
# LC_ALL=C so that sort's byte order and awk's string comparison agree -- gawk uses strcoll
# for `<` outside the C locale, and a collation disagreeing with `sort` would misalign the
# streams. IDs are `<sample_set>.<caller>_<n>`, never numeric, so there is no
# lexical-vs-numeric trap; the "" concatenations below keep every comparison in string
# context regardless.
#
# Getting this wrong does not error, it UNDER-JOINS: depth goes silently missing for some
# records and the downstream depth filters then run on padding. That is precisely the class of
# failure this file exists to prevent, so both streams are checked for ascending ID.

# The depth table is read with getline rather than as a main input, so the VCF can be
# streamed against it; blanking ARGV[1] removes it from awk's main input list. It is still
# identified by ARGV position, not by the FNR==NR idiom: when no sample in the cohort has an
# alignment the depth table is legitimately EMPTY, and FNR==NR then never becomes false --
# awk reads the VCF as if it were the depth table and emits nothing at all. A cohort with no
# alignments is a supported configuration, not an error, so it has to come out the other
# side as a valid VCF with padded columns.
BEGIN {
    FS = OFS = "\t"
    depth_file = ARGV[1]
    ARGV[1] = ""
    have = next_depth()
}

# Pulls the next usable depth row into d_id / d_sample / d_vals. Returns 0 at end of stream.
function next_depth(   r, line, D) {
    while (1) {
        r = (getline line < depth_file)
        # -1, not 0: an unreadable table would otherwise look like an empty one and pad every
        # column silently, which is indistinguishable from a cohort with no alignments.
        if (r < 0) {
            printf("ERROR: cannot read depth table '%s'\n", depth_file) > "/dev/stderr"
            exit 1
        }
        if (r == 0) return 0
        split(line, D, "\t")
        if (D[1] == "ID" || D[1] ~ /^#/) continue      # tolerate a header line
        if (D[1] "" < d_prev "") {
            printf("ERROR: depth table is not sorted by ID (%s follows %s); " \
                   "merge_depth.awk needs both inputs ID-sorted\n",
                   D[1], d_prev) > "/dev/stderr"
            exit 1
        }
        d_prev = D[1] ""
        d_id = D[1] ""
        d_sample = D[2] ""
        d_vals = D[3] ":" D[4] ":" D[5] ":" D[6]
        return 1
    }
    return 0
}

# --- header ------------------------------------------------------------------------
/^##/ { print; next }

/^#CHROM/ {
    print "##FORMAT=<ID=DHFFC,Number=1,Type=Float,Description=\"duphold depth fold-change vs flanks\">"
    print "##FORMAT=<ID=DHBFC,Number=1,Type=Float,Description=\"duphold depth fold-change vs GC-matched bins\">"
    # Per duphold: fold-change vs the rest of the CHROMOSOME the variant is on, not the genome.
    # The distinction matters -- a chromosome-relative value cannot reveal whole-chromosome
    # ploidy, because the denominator is halved on a hemizygous chromosome too.
    print "##FORMAT=<ID=DHFC,Number=1,Type=Float,Description=\"duphold depth fold-change vs the median depth of the rest of the chromosome\">"
    print "##FORMAT=<ID=DHBZ,Number=1,Type=Float,Description=\"duphold depth z-score vs GC-matched bins\">"
    for (i = 10; i <= NF; i++) sample_name[i] = $i
    print
    next
}

# --- records ----------------------------------------------------------------------
{
    id = $3 ""
    if (seen && id < v_prev "") {
        printf("ERROR: cohort VCF body is not sorted by ID (%s follows %s); " \
               "merge_depth.awk needs both inputs ID-sorted\n", id, v_prev) > "/dev/stderr"
        exit 1
    }
    v_prev = id
    seen = 1

    # Rebuilt only when the ID changes, so two adjacent records sharing an ID both get the
    # same block rather than the second one silently getting nothing.
    if (!built || id != cur_id) {
        delete cur
        while (have && d_id < id) have = next_depth()   # depth for a record not in the VCF
        while (have && d_id == id) {
            cur[d_sample] = d_vals
            have = next_depth()
        }
        cur_id = id
        built = 1
    }

    $9 = $9 ":DHFFC:DHBFC:DHFC:DHBZ"
    for (i = 10; i <= NF; i++)
        $i = $i ":" ((sample_name[i] in cur) ? cur[sample_name[i]] : ".:.:.:.")
    print
}
