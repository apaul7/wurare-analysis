# Rare / potentially-damaging subsets of the AnnotSV TSV.
#
# Usage:
#   awk -v TIER=1 -v RARE_AF=0.01 -v ACMG_MIN=4 -v PLI_MIN=0.9 -v LOEUF_MAX=2 \
#       -v RANK_MIN=0.9 -f filter_annotsv_tsv.awk annotsv.tsv > tier1.tsv
#
# This is the one place in the pipeline that REMOVES rows, and it is allowed to only because
# the unfiltered TSV is published beside it and never replaced. The "annotate hard, filter
# softly" rule is about the callset; a report is a view of it, not the thing itself.
#
# BOTH TIERS KEEP SPLIT ROWS ONLY. AnnotSV runs in `-annotationMode both`, so the unfiltered
# TSV carries full rows (one per SV) and split rows (one per SV x gene). The gene-level
# columns every criterion below depends on -- Location, Overlapped_CDS_percent, GnomAD_pLI,
# LOEUF_bin, OMIM_morbid -- are per-transcript, and AnnotSV leaves them BLANK on full rows.
# Filtering a full row on coding-ness is therefore not strict but inert, which is the whole
# reason the mode changed. The consequence to keep in mind when reading the output: a row
# here is an SV x GENE, so one SV spanning four genes contributes up to four rows, and the
# row count is not a variant count. Deduplicate on the ID column for that.
#
# Two tiers, because one threshold cannot serve both questions a reviewer asks:
#
#   TIER 1 -- "is there a single answer in this cohort". PASS-only, coding, ACMG 4/5, in a
#   constrained gene, not seen in AnnotSV's benign SV sets. EXPECT ZERO TO A HANDFUL OF ROWS
#   PER SAMPLE, AND OFTEN NONE AT ALL. An empty tier-1 file is the normal result, not a bug.
#
#   TIER 2 -- the working list. Drops only what population frequency says is common, keeps
#   anything coding or in an OMIM morbid gene, and asks for ACMG 3+ or a high ranking score.
#   Hundreds to thousands of rows. NO_CALLER_SUPPORT and DEPTH_UNSUPPORTED records are kept
#   deliberately, with the tag still visible in the FILTER column -- a caller missing an event
#   is weak evidence of absence.
#
# Every column is located BY NAME from the header line, never by position: AnnotSV's column
# set shifts between releases and with the annotation bundle, and a filter reading the wrong
# column is worse than no filter at all. A column a tier needs but cannot find is named in the
# header comment line rather than guessed at.
#
# One caveat that belongs here rather than in a review meeting. Tier 1's reliance on
# FILTER == PASS means DEPTH_UNSUPPORTED only excludes anything once duphold has actually run
# -- which needs alignments, and on chrX/chrY also needs --somalier_sites for the ploidy the
# threshold is scaled by. Without those, that criterion silently passes everything.

function present(v) { return (v != "" && v != "." && v != "NA") }

function col(name) { return (name in c) ? c[name] : 0 }

function get(name,   i) { i = col(name); return (i == 0) ? "" : $i }

# AnnotSV writes several gene-based columns as a list when one row spans many genes. Take the
# worst case: the highest pLI, the lowest (most constrained) LOEUF bin, the highest benign AF.
# A scalar field is a one-element list, so this covers both shapes.
function extreme(s, want_max,   n, i, parts, best, have, v) {
    n = split(s, parts, /[;,|]/)
    have = 0
    for (i = 1; i <= n; i++) {
        if (parts[i] !~ /^-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$/) continue
        v = parts[i] + 0
        if (!have || (want_max ? v > best : v < best)) { best = v; have = 1 }
    }
    return have ? best : ""
}

# "4", "full=4" and "NA" all occur depending on release and row. Pull the digit rather than
# comparing the whole cell.
function acmg(   s) {
    s = get("ACMG_class")
    return match(s, /[1-5]/) ? substr(s, RSTART, 1) : ""
}

function has_tag(tag) { return get("FILTER") ~ ("(^|;)" tag "($|;)") }

# Split rows only reach this, so both columns are populated and a blank is a real negative
# rather than the "not applicable" it means on a full row.
function is_coding(   pct) {
    if (get("Location") ~ /CDS/) return 1
    pct = extreme(get("Overlapped_CDS_percent"), 1)
    return (pct != "" && pct > 0)
}

# The row-type gate, applied to both tiers. AnnotSV runs in `both` mode, so a full row is the
# same SV summarised without its gene-level annotation -- keeping it would duplicate every
# variant AND carry blanks through every criterion below. If the column is missing entirely
# the gate cannot be applied; that is reported in the header line rather than guessed at.
function is_split() {
    return (col("Annotation_mode") == 0 || get("Annotation_mode") == "split")
}

function rule_text() {
    if (TIER == 1)
        return "tier1 strict, one row per SV x gene: Annotation_mode=split; FILTER=PASS;" \
               " coding (Location~CDS or Overlapped_CDS_percent>0);" \
               " ACMG_class>=" ACMG_MIN "; GnomAD_pLI>=" PLI_MIN " or LOEUF_bin<=" LOEUF_MAX ";" \
               " B_gain_AFmax and B_loss_AFmax <" RARE_AF " or blank"
    return "tier2 loose, one row per SV x gene: Annotation_mode=split; no COMMON_GNOMAD or" \
           " COMMON_INTERNAL tag (other FILTER tags kept and visible); coding or OMIM_morbid;" \
           " ACMG_class>=" ACMG_MIN " or AnnotSV_ranking_score>=" RANK_MIN
}

BEGIN {
    FS = OFS = "\t"
    if (TIER == "")      TIER = 1
    if (RARE_AF == "")   RARE_AF = 0.01
    if (ACMG_MIN == "")  ACMG_MIN = (TIER == 1 ? 4 : 3)
    if (PLI_MIN == "")   PLI_MIN = 0.9
    if (LOEUF_MAX == "") LOEUF_MAX = 2
    if (RANK_MIN == "")  RANK_MIN = 0.9
}

NR == 1 {
    for (i = 1; i <= NF; i++) c[$i] = i

    # Which columns each tier actually reads. Named here so a release that renames one shows
    # up as a line in the file, rather than as a criterion that quietly stopped applying.
    need = (TIER == 1) \
        ? "FILTER Annotation_mode Location Overlapped_CDS_percent ACMG_class GnomAD_pLI LOEUF_bin B_gain_AFmax B_loss_AFmax" \
        : "FILTER Annotation_mode Location Overlapped_CDS_percent OMIM_morbid ACMG_class AnnotSV_ranking_score"
    n = split(need, want, " ")
    miss = ""
    for (i = 1; i <= n; i++)
        if (!(want[i] in c)) miss = miss (miss == "" ? "" : ",") want[i]

    print "## " rule_text() (miss == "" ? "" : "  | MISSING COLUMNS, criterion not applied: " miss)
    print
    next
}

TIER == 1 {
    if (!is_split()) next
    if (get("FILTER") != "PASS") next
    if (!is_coding()) next

    a = acmg()
    if (a == "" || a + 0 < ACMG_MIN + 0) next

    pli = extreme(get("GnomAD_pLI"), 1)
    loeuf = extreme(get("LOEUF_bin"), 0)
    if (!((pli != "" && pli >= PLI_MIN + 0) || (loeuf != "" && loeuf <= LOEUF_MAX + 0))) next

    # Blank means AnnotSV found no benign SV overlapping this one, which IS the rare case --
    # so blank passes, and only a measured AF at or above the threshold drops the row.
    bgain = extreme(get("B_gain_AFmax"), 1)
    bloss = extreme(get("B_loss_AFmax"), 1)
    if (bgain != "" && bgain >= RARE_AF + 0) next
    if (bloss != "" && bloss >= RARE_AF + 0) next

    print
    next
}

{
    if (!is_split()) next
    if (has_tag("COMMON_GNOMAD") || has_tag("COMMON_INTERNAL")) next
    # Coding OR in a disease gene, so a promoter or regulatory hit in an OMIM morbid gene
    # survives -- that is the whole point of the second tier.
    if (!(is_coding() || get("OMIM_morbid") ~ /yes/)) next

    a = acmg()
    rank = extreme(get("AnnotSV_ranking_score"), 1)
    if (!((a != "" && a + 0 >= ACMG_MIN + 0) || (rank != "" && rank >= RANK_MIN + 0))) next

    print
}
