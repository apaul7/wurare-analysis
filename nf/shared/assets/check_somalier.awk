# Per-sample ploidy and the PED cross-check, from somalier's own output.
#
# Usage:
#   awk -v SAMPLES=cohort.samples.tsv -v PAIRS=cohort.pairs.tsv -v PED=cohort.ped \
#       -v OUT_PLOIDY=ploidy.tsv -v OUT_PED=cohort.inferred.ped \
#       -f check_somalier.awk </dev/null
#
# Writes:
#   ploidy.tsv           sample <tab> x_copies <tab> y_copies <tab> ped_sex <tab> agreement
#   cohort.inferred.ped  the input PED with column 5 replaced by the inferred sex
#   warnings to stderr
#
# Everything happens in BEGIN and there is no main input: three files are read by name, and
# two are written. Reading them positionally with FNR==NR chains does not survive a third
# file, and there is no stream to run against anyway.
#
# --- why ploidy is derived from depth, not read off somalier's `sex` column -----------
# `somalier relate --infer` can write an inferred sex into samples.tsv, but which column
# carries it -- and whether it holds the pedigree value or the inferred one -- depends on how
# the run was invoked. The two depth columns do not: X_depth_mean and Y_depth_mean are
# comparable against depth_mean, which is dominated by autosomal sites. One X against two
# reads ~0.5 there, and that ratio is what the somalier HTML plots sex on. Deriving it here
# means this file depends on three stable column NAMES rather than on an inference flag's
# output shape -- and the names are asserted, not assumed, because a renamed column would
# otherwise read as a cohort of unknown-sex samples and silently restore the old
# autosome-only behaviour.
#
# --- copies, not sex ------------------------------------------------------------------
# x_copies and y_copies are emitted separately and are what tag_filters.awk consumes. A
# karyotype that is neither XX nor XY -- XXY, X0 -- still gets usable per-chromosome copy
# numbers even though its `agreement` is UNKNOWN, because the depth threshold on chrX only
# needs to know how many X copies to expect. "." means undetermined: an ambiguous ratio is
# left undetermined rather than rounded to the nearer karyotype, and tag_filters then skips
# that sample exactly as it skips a sample with no depth measured at all.
#
# --- the PED is the cross-check, never the source -------------------------------------
# A disagreement is reported and the DATA wins: sample swaps and mislabelled PEDs are why
# this stage exists. It is deliberately not fatal -- an ambiguous ratio is a property of some
# real samples, not an error -- so it is made loud instead: here on stderr, in ploidy.tsv's
# agreement column, and in a ##SAMPLE_SEX header line that tag_filters.awk writes into the
# published VCF.

function fail(msg) {
    printf "check_somalier: ERROR %s\n", msg > "/dev/stderr"
    exit 1
}

# Column indices by NAME. A somalier release that renames or reorders a column must fail
# here, naming it, rather than reading a number out of the wrong field.
function header_index(line, want,    n, i, parts, name) {
    n = split(line, parts, "\t")
    for (i = 1; i <= n; i++) {
        name = parts[i]
        sub(/^#/, "", name)
        if (name == want) return i
    }
    return 0
}

function ratio(num, den) { return (den + 0 > 0) ? (num + 0) / (den + 0) : -1 }

BEGIN {
    FS = OFS = "\t"

    if (SAMPLES == "")    fail("SAMPLES is required (somalier relate's *.samples.tsv)")
    if (PED == "")        fail("PED is required")
    if (OUT_PLOIDY == "") fail("OUT_PLOIDY is required")
    if (OUT_PED == "")    fail("OUT_PED is required")

    # Cut points, overridable with -v. Defaults sit either side of the two real modes (~0.5
    # and ~1.0 for X; ~0 and ~0.4 for Y -- chrY is repeat-rich and much of it unmappable, so
    # a male never reaches 1.0 there). The gap between each pair is deliberately wide: what
    # falls in it is undetermined, and undetermined costs a filter rather than corrupting one.
    if (X_HEMI_MAX == "")    X_HEMI_MAX = 0.65
    if (X_DIPLOID_MIN == "") X_DIPLOID_MIN = 0.80
    if (Y_ABSENT_MAX == "")  Y_ABSENT_MAX = 0.05
    if (Y_PRESENT_MIN == "") Y_PRESENT_MIN = 0.15

    read_samples()
    read_ped()
    write_ploidy()
    write_inferred_ped()
    if (PAIRS != "") check_pairs()
}

function read_samples(    line, hdr, i_sample, i_depth, i_x, i_y, n, parts, s, xr, yr) {
    if ((getline line < SAMPLES) <= 0)
        fail("cannot read " SAMPLES " -- somalier relate produced no samples table")
    hdr = line

    i_sample = header_index(hdr, "sample_id")
    i_depth  = header_index(hdr, "depth_mean")
    i_x      = header_index(hdr, "X_depth_mean")
    i_y      = header_index(hdr, "Y_depth_mean")
    if (!i_sample) fail("no sample_id column in " SAMPLES)
    if (!i_depth)  fail("no depth_mean column in " SAMPLES)
    if (!i_x)      fail("no X_depth_mean column in " SAMPLES)
    if (!i_y)      fail("no Y_depth_mean column in " SAMPLES)

    n_seen = 0
    while ((getline line < SAMPLES) > 0) {
        if (line == "" || line ~ /^#/) continue
        n = split(line, parts, "\t")
        if (n < i_y) fail("short row in " SAMPLES ": " line)
        s = parts[i_sample]

        # Zero depth means somalier genotyped nothing for this sample. The usual cause is a
        # sites file built against a different reference than the alignments, which otherwise
        # produces a complete and entirely wrong table. Fatal: no ploidy call made from it
        # would mean anything, and the run must not proceed quietly to the depth filter.
        if (parts[i_depth] + 0 <= 0)
            fail("sample " s " has depth_mean " parts[i_depth] " -- somalier genotyped no " \
                 "sites for it. The likeliest cause is a --somalier_sites file built " \
                 "against a different reference than --alignment_reference")

        xr = ratio(parts[i_x], parts[i_depth])
        yr = ratio(parts[i_y], parts[i_depth])

        x_copies[s] = (xr < 0) ? "." : \
                      (xr <= X_HEMI_MAX + 0) ? 1 : \
                      (xr >= X_DIPLOID_MIN + 0) ? 2 : "."
        y_copies[s] = (yr < 0) ? "." : \
                      (yr <= Y_ABSENT_MAX + 0) ? 0 : \
                      (yr >= Y_PRESENT_MIN + 0) ? 1 : "."
        x_ratio[s] = xr
        y_ratio[s] = yr
        seen[s] = 1
        n_seen++
        order[n_seen] = s
    }
    close(SAMPLES)
    if (n_seen == 0) fail(SAMPLES " has a header but no sample rows")
}

# The PED is read in file order and both outputs follow it, so they are byte-stable across
# runs -- an unstable QC table would re-hash every task that reads it and defeat -resume, the
# trap the sorted collects elsewhere in this pipeline exist for.
function read_ped(    line, n, parts) {
    n_ped = 0
    while ((getline line < PED) > 0) {
        n_ped++
        ped_line[n_ped] = line
        if (line ~ /^[ \t]*#/ || line ~ /^[ \t]*$/) continue
        n = split(line, parts, "[ \t]+")
        if (n < 6) fail("PED line " n_ped " has " n " columns, expected 6: " line)
        ped_sample[n_ped] = parts[2]
        ped_sex[parts[2]] = parts[5]
    }
    close(PED)
}

# male=1, female=2, 0 undetermined. Only the two complete karyotypes get a sex: an XXY or X0
# sample has usable per-chromosome copy numbers but no answer to this question, and inventing
# one would write a wrong value into the PED that sex-stratified AF then uses.
function inferred_sex(s) {
    if (x_copies[s] == 1 && y_copies[s] == 1) return 1
    if (x_copies[s] == 2 && y_copies[s] == 0) return 2
    return 0
}

function agreement(s,    inf, ped) {
    inf = inferred_sex(s)
    ped = (s in ped_sex) ? ped_sex[s] + 0 : -1
    if (inf == 0) return "UNKNOWN"
    if (ped != 1 && ped != 2) return "PED_MISSING"
    return (inf == ped) ? "AGREES" : "DISAGREES"
}

function write_ploidy(    i, s, agr, ped) {
    printf "#sample\tx_copies\ty_copies\tped_sex\tagreement\n" > OUT_PLOIDY
    for (i = 1; i <= n_seen; i++) {
        s = order[i]
        agr = agreement(s)
        ped = (s in ped_sex) ? ped_sex[s] : "."
        printf "%s\t%s\t%s\t%s\t%s\n", s, x_copies[s], y_copies[s], ped, agr > OUT_PLOIDY

        if (agr == "DISAGREES")
            printf "check_somalier: WARNING sample %s -- PED says sex=%s, the alignments say " \
                   "%s (X depth ratio %.2f, Y %.2f). The DATA is used; check for a sample " \
                   "swap or a mislabelled PED row\n", \
                   s, ped, (inferred_sex(s) == 1 ? "male" : "female"), \
                   x_ratio[s], y_ratio[s] > "/dev/stderr"
        else if (agr == "UNKNOWN")
            printf "check_somalier: NOTE sample %s has no determined karyotype " \
                   "(X depth ratio %.2f, Y %.2f); chrX/chrY depth tagging is skipped for it\n", \
                   s, x_ratio[s], y_ratio[s] > "/dev/stderr"
        else if (agr == "PED_MISSING")
            printf "check_somalier: NOTE sample %s has no usable PED sex; the inferred one " \
                   "is used\n", s > "/dev/stderr"
    }
    close(OUT_PLOIDY)
}

# The input PED with column 5 replaced where a sex was determined. Family and parent columns
# are never touched: somalier is asked what the data says about sex, not to invent pedigrees.
# A sample with no determined sex keeps whatever the PED had, including 0.
function write_inferred_ped(    i, s, inf, n, parts, j, out) {
    for (i = 1; i <= n_ped; i++) {
        s = ped_sample[i]
        inf = (s != "" && (s in seen)) ? inferred_sex(s) : 0
        if (inf == 0) {
            print ped_line[i] > OUT_PED
            continue
        }
        n = split(ped_line[i], parts, "[ \t]+")
        parts[5] = inf
        out = parts[1]
        for (j = 2; j <= n; j++) out = out "\t" parts[j]
        print out > OUT_PED
    }
    close(OUT_PED)
}

# Relatedness against what the PED implies. Warning only, and deliberately so: this is a
# sample-swap detector for the Talos inheritance tail, and a cohort is not wrong to contain a
# pair whose true relationship nobody recorded. A pair somalier gives no expectation for is
# skipped rather than guessed at.
function check_pairs(    line, hdr, i_a, i_b, i_rel, i_exp, n, parts, rel, want) {
    if ((getline line < PAIRS) <= 0) {
        printf "check_somalier: NOTE %s is empty; relatedness not checked\n", \
               PAIRS > "/dev/stderr"
        return
    }
    hdr = line
    i_a   = header_index(hdr, "sample_a")
    i_b   = header_index(hdr, "sample_b")
    i_rel = header_index(hdr, "relatedness")
    i_exp = header_index(hdr, "expected_relatedness")
    if (!i_a || !i_b || !i_rel) fail("no sample_a/sample_b/relatedness columns in " PAIRS)
    if (!i_exp) {
        printf "check_somalier: NOTE %s has no expected_relatedness column, so pedigree " \
               "relatedness is not checked\n", PAIRS > "/dev/stderr"
        close(PAIRS)
        return
    }

    while ((getline line < PAIRS) > 0) {
        if (line == "" || line ~ /^#/) continue
        n = split(line, parts, "\t")
        if (n < i_exp) continue
        want = parts[i_exp] + 0
        rel = parts[i_rel] + 0
        if (want < 0) continue          # somalier writes -1 where the PED implies nothing
        if (want >= 0.45 && rel < 0.35)
            printf "check_somalier: WARNING %s and %s are declared related (expected %.2f) " \
                   "but measure %.2f -- possible sample swap or wrong PED\n", \
                   parts[i_a], parts[i_b], want, rel > "/dev/stderr"
        else if (want < 0.10 && rel > 0.20)
            printf "check_somalier: WARNING %s and %s are declared unrelated (expected %.2f) " \
                   "but measure %.2f -- possible duplicate sample or unrecorded relationship\n", \
                   parts[i_a], parts[i_b], want, rel > "/dev/stderr"
    }
    close(PAIRS)
}
