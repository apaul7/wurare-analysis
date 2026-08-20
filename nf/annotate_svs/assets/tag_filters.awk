# Soft filter tags. Annotate hard, filter softly.
#
# Usage:
#   awk -v POP_AF=0.01 -v INT_AF=0.03 -v DEL_DHFFC=0.7 -v DUP_DHBFC=1.3 -v MIN_CALLERS=2 \
#       -v PLOIDY=ploidy.tsv -v BUILD=GRCh38 \
#       -f tag_filters.awk cohort.vcf
#
# Every record is KEPT. Failing a criterion writes a tag into the FILTER column and nothing
# else -- no record is ever removed. One run then serves Talos, manual review and a QC
# report without re-running the expensive half, and each consumer sets its own threshold
# from the tags rather than from a filtered file it cannot widen again.
#
# Tags written:
#   COMMON_GNOMAD       population AF at or above POP_AF in either gnomAD callset
#   COMMON_INTERNAL     internal cohort AF at or above INT_AF
#   DEPTH_UNSUPPORTED   a DEL whose depth did not drop, or a DUP whose depth did not rise
#   NO_CALLER_SUPPORT   fewer than MIN_CALLERS callers found this record
#
# --- the sex chromosomes -------------------------------------------------------------
# The DUP threshold is scaled by the number of copies the SAMPLE is expected to carry at the
# LOCUS, which is why this file takes a ploidy table. duphold's DHBFC is normalized against
# GC-matched bins GENOME-wide, so a hemizygous male chrX sits near 0.5 and a real 2-copy
# duplication there reads about 1.0 -- under the flat 1.3 threshold. Judged against
# 1.3 x 1/2 = 0.65 it reads as supported, which is what it is.
#
# The DEL threshold is NOT scaled. DHFFC is normalized against the variant's own FLANKS,
# which are equally hemizygous, so a male chrX region reads about 1.0 at baseline and about 0
# when deleted; the 0.7 threshold already means the right thing there. (This file said the
# opposite for a long time and was wrong.)
#
# Ploidy comes from the alignments (somalier, via nf/shared/assets/check_somalier.awk), never from the
# PED -- PED sex cannot be assumed correct. Where the two disagree, the data is
# used and the disagreement is written into a ##SAMPLE_SEX header line below, so it reaches
# every consumer of this VCF rather than living in a log nobody reads.
#
# With no PLOIDY table, or for a sample missing from it, or one whose karyotype came back
# undetermined, the sex chromosomes stay exempt exactly as they were before any of this
# existed. Absence of evidence is not evidence.
#
# PAR is diploid in every karyotype, so a record overlapping one gets the unscaled thresholds
# for everyone. Overlap rather than containment, deliberately: a record straddling a PAR
# boundary is judged as diploid, which is the LENIENT direction -- a hemizygous baseline
# tested against the diploid threshold under-tags, where the reverse would tag true calls.
# The coordinates are GRCh38's, this pipeline's only build; any other value of BUILD falls
# back to exempting the sex chromosomes entirely and says so, because PAR coordinates from
# the wrong build would be silent.
#
# BND and INS are exempt too: a breakend has no interval to measure depth over, and an
# insertion has no reference span to lose or gain coverage in. duphold still emits numbers
# for them; they just do not mean what the threshold assumes.
#
# A record with no depth measured at all -- every sample lacking an alignment -- is not
# tagged. Absence of evidence is not evidence, and the design is explicit that depth and caller
# corroboration are score inputs rather than gates.

function info_value(info, key,    n, i, parts) {
    n = split(info, parts, ";")
    for (i = 1; i <= n; i++) {
        if (parts[i] ~ ("^" key "=")) return substr(parts[i], length(key) + 2)
    }
    return ""
}

function present(v) { return (v != "" && v != ".") }

# INFO fields declared Number=A carry one value per ALT allele, so a multi-allelic record
# holds a comma list. Bare "+0" takes the leading number, which silently evaluates the record
# on its first allele alone -- a record whose SECOND allele is common at 0.9 reads as 0.001
# and escapes the tag. Multi-allelic records are a live path here: symbolic_alt.awk
# deliberately preserves them. Worst case across alleles is the honest reading of "is this
# common".
function max_value(v,    n, i, parts, best, have) {
    n = split(v, parts, ",")
    have = 0
    for (i = 1; i <= n; i++) {
        if (!present(parts[i])) continue
        if (!have || parts[i] + 0 > best) { best = parts[i] + 0; have = 1 }
    }
    return have ? best : ""
}

# Internal cohort AF, counted over EVERY sample rather than over called alleles only.
#
# Deliberately not the AF bcftools +fill-tags wrote. fill-tags computes AN over called alleles
# and SVDB leaves a non-calling sample as "./." rather than "0/0" -- measured, asserted in
# tests/test_svdb_assumptions.py -- so a private variant scored AF=0.5 at any cohort size and
# COMMON_INTERNAL fired on almost the whole callset.
#
# ponytail: a no-call counts as reference, understating AF for a shallow sample. Excluding
# them, as fill-tags does, overstates it far worse at this cohort size. Revisit when the
# frozen cross-run frequency database lands.
function cohort_af(gt_idx,    c, n_alt, n_alleles, g, parts, np, j) {
    n_alt = 0; n_alleles = 0
    for (c = 10; c <= NF; c++) {
        split($c, parts, ":")
        g = parts[gt_idx]
        gsub(/\|/, "/", g)
        np = split(g, alleles, "/")
        for (j = 1; j <= np; j++) {
            # "." is a no-call for that chromosome copy; it still occupies a slot in the
            # denominator, which is the whole point of counting over all samples.
            n_alleles += 1
            if (alleles[j] != "." && alleles[j] + 0 > 0) n_alt += 1
        }
    }
    return (n_alleles > 0) ? n_alt / n_alleles : ""
}

# PAR coordinates, 1-based inclusive, from the build's own documentation. Held as parallel
# arrays rather than a BED asset: there are four intervals, they are constants of the
# reference and not of a run, and a file would be one more thing to stage and to get wrong.
#
# GRCh38 only, because this pipeline is GRCh38 only. Any other build gets no coordinates and
# the sex chromosomes stay exempt -- coordinates from the wrong build would be silent, and
# a table nobody runs against is a table nobody notices is wrong.
function load_par(    b) {
    b = toupper(BUILD)
    if (b == "GRCH38" || b == "HG38") {
        par_n = 4
        par_chr[1] = "X"; par_beg[1] = 10001;     par_end[1] = 2781479
        par_chr[2] = "X"; par_beg[2] = 155701383; par_end[2] = 156030895
        par_chr[3] = "Y"; par_beg[3] = 10001;     par_end[3] = 2781479
        par_chr[4] = "Y"; par_beg[4] = 56887903;  par_end[4] = 57217415
    }
    else {
        par_n = 0
    }
}

# sample -> expected copies on X and on Y. "." means undetermined, and an undetermined sample
# is skipped rather than assumed diploid: assuming would silently apply the autosomal DUP
# threshold to a hemizygous sample, which is the defect this whole file is undoing.
function load_ploidy(    line, n, parts) {
    have_ploidy = 0
    if (PLOIDY == "") return
    while ((getline line < PLOIDY) > 0) {
        if (line ~ /^#/ || line ~ /^[ \t]*$/) continue
        n = split(line, parts, "\t")
        if (n < 3) continue
        x_copies[parts[1]] = parts[2]
        y_copies[parts[1]] = parts[3]
        if (n >= 5) sample_agreement[parts[1]] = parts[5]
        if (n >= 4) sample_ped_sex[parts[1]] = parts[4]
        have_ploidy = 1
    }
    close(PLOIDY)
}

BEGIN {
    FS = OFS = "\t"
    if (POP_AF == "")      POP_AF = 0.01
    if (INT_AF == "")      INT_AF = 0.03
    if (DEL_DHFFC == "")   DEL_DHFFC = 0.7
    if (DUP_DHBFC == "")   DUP_DHBFC = 1.3
    if (MIN_CALLERS == "") MIN_CALLERS = 2

    load_ploidy()
    load_par()
    if (have_ploidy && par_n == 0)
        printf "tag_filters: NOTE no PAR coordinates for genome build '%s', so chrX/chrY " \
               "stay exempt from the depth tag. This pipeline is GRCh38 only\n", \
               BUILD > "/dev/stderr"
}

# Bare chromosome name: chr-prefixed and bare references both reach this pipeline, and the
# PAR table is written once.
function bare(chrom,    c) { c = chrom; sub(/^chr/, "", c); return c }

function in_par(chrom, beg, end,    c, i) {
    c = bare(chrom)
    for (i = 1; i <= par_n; i++)
        if (par_chr[i] == c && beg <= par_end[i] && end >= par_beg[i]) return 1
    return 0
}

# Expected copies for one sample over one record: 2 on an autosome or in PAR, the sample's own
# X or Y count outside it, -1 where nothing is known and the sample must be skipped.
function expected_copies(chrom, beg, end, sample,    c, v) {
    c = bare(chrom)
    if (c != "X" && c != "Y") return 2
    if (!have_ploidy || par_n == 0) return -1
    if (in_par(chrom, beg, end)) return 2
    if (!(sample in x_copies)) return -1
    v = (c == "X") ? x_copies[sample] : y_copies[sample]
    return (v == "." || v == "") ? -1 : v + 0
}

/^#CHROM/ {
    n_samples = NF - 9
    for (c = 10; c <= NF; c++) sample_name[c] = $c

    # Whether COMMON_INTERNAL can mean anything is a property of the cohort, decided once
    # here. Below the size where the floor 1/(2N) already clears the threshold, every carrier
    # is tagged and the tag means "is a carrier" rather than "is common" -- worse than no tag,
    # because a reviewer reads it as a reason to deprioritize. So it is not written, and the
    # run says so. At ~17 samples (11 trios) 1/34 = 0.029 sits just under the 0.03 default, so
    # a singleton passes and a doubleton is tagged, which is the intent.
    int_af_usable = (n_samples > 0 && 1.0 / (2 * n_samples) < INT_AF + 0)
    if (!int_af_usable) {
        printf "tag_filters: NOTE %d sample(s) cannot express a frequency below %s " \
               "(floor is 1/(2N) = %.4g), so COMMON_INTERNAL is not written for this run\n", \
               n_samples, INT_AF, (n_samples > 0 ? 1.0 / (2 * n_samples) : 1) > "/dev/stderr"
    }

    print "##INFO=<ID=INTERNAL_AF,Number=1,Type=Float,Description=\"Internal cohort allele " \
          "frequency over all samples, counting no-calls as reference. Not INFO/AF, which " \
          "bcftools +fill-tags computes over called alleles only\">"
    print "##FILTER=<ID=COMMON_GNOMAD,Description=\"Population AF >= " POP_AF " in gnomAD SV or CNV\">"
    print "##FILTER=<ID=COMMON_INTERNAL,Description=\"Internal cohort AF >= " INT_AF "\">"
    print "##FILTER=<ID=DEPTH_UNSUPPORTED,Description=\"duphold depth does not support the call type; DEL/DUP/CNV only, with the DUP threshold scaled by the sample's expected copy number on chrX/chrY\">"
    print "##FILTER=<ID=NO_CALLER_SUPPORT,Description=\"Fewer than " MIN_CALLERS " callers support this record\">"

    # The ploidy each sample was filtered under, in the file it filtered. A reviewer asking
    # "why was this male chrX duplication kept" -- or "why does this sample's sex not match
    # the referral" -- can answer it from the VCF alone, without the QC directory. agreement
    # is PED vs data: DISAGREES means the PED said otherwise and the data was used.
    for (c = 10; c <= NF; c++) {
        if (!(sample_name[c] in x_copies)) continue
        printf "##SAMPLE_SEX=<ID=%s,x_copies=%s,y_copies=%s,ped_sex=%s,agreement=%s>\n", \
               sample_name[c], x_copies[sample_name[c]], y_copies[sample_name[c]], \
               (sample_name[c] in sample_ped_sex) ? sample_ped_sex[sample_name[c]] : ".", \
               (sample_name[c] in sample_agreement) ? sample_agreement[sample_name[c]] : "."
    }
    print
    next
}

/^#/ { print; next }

{
    tags = ""

    # --- population frequency ---------------------------------------------------
    sv_af = max_value(info_value($8, "gnomad_sv_AF"))
    cnv_af = max_value(info_value($8, "gnomad_cnv_AF"))
    if ((present(sv_af)  && sv_af  + 0 >= POP_AF + 0) ||
        (present(cnv_af) && cnv_af + 0 >= POP_AF + 0)) {
        tags = "COMMON_GNOMAD"
    }

    # --- internal cohort frequency ----------------------------------------------
    # From the genotypes, not INFO/AF -- see cohort_af above. Still a floor at small N, which
    # is why the design leans on the population AF while the cohort is small. Written back as
    # INTERNAL_AF; INFO/AF is left alone because Talos reads it.
    gt_idx = 0
    n = split($9, keys, ":")
    for (i = 1; i <= n; i++) if (keys[i] == "GT") gt_idx = i

    int_af = (gt_idx > 0 && NF >= 10) ? cohort_af(gt_idx) : ""
    if (present(int_af)) {
        # Written even when the cohort is too small for the tag: the number is still true.
        # What small N cannot support is the judgement that it is high.
        $8 = $8 ";INTERNAL_AF=" sprintf("%.6g", int_af)
        if (int_af_usable && int_af + 0 >= INT_AF + 0) {
            tags = (tags == "") ? "COMMON_INTERNAL" : tags ";COMMON_INTERNAL"
        }
    }

    # --- caller corroboration ---------------------------------------------------
    # Flag, never drop: a caller missing an event is weak evidence of absence, because the
    # callers have different sensitivity profiles.
    ncaller = info_value($8, "NCALLER")
    if (present(ncaller) && ncaller + 0 < MIN_CALLERS + 0) {
        tags = (tags == "") ? "NO_CALLER_SUPPORT" : tags ";NO_CALLER_SUPPORT"
    }

    # --- depth -------------------------------------------------------------------
    svtype = info_value($8, "SVTYPE")
    if (svtype == "DEL" || svtype == "DUP" || svtype == "CNV") {
        n = split($9, keys, ":")
        idx_ffc = 0; idx_bfc = 0
        for (i = 1; i <= n; i++) {
            if (keys[i] == "DHFFC") idx_ffc = i
            if (keys[i] == "DHBFC") idx_bfc = i
        }

        # END for the PAR overlap test. A record whose END is absent or nonsensical is
        # measured over POS alone -- normalize_records.awk has already rejected END < POS.
        rec_end = info_value($8, "END")
        if (!present(rec_end) || rec_end + 0 < $2 + 0) rec_end = $2

        supported = 0
        measured = 0
        for (c = 10; c <= NF; c++) {
            copies = expected_copies($1, $2 + 0, rec_end + 0, sample_name[c])
            # -1 is "nothing known here" and 0 is "no copies to gain or lose" -- a chrY call
            # in a sample with no Y. Neither yields a threshold that means anything, so the
            # sample contributes nothing rather than a negative opinion.
            if (copies <= 0) continue
            split($c, vals, ":")
            if (svtype == "DEL") {
                if (idx_ffc > 0 && present(vals[idx_ffc])) {
                    measured = 1
                    # Unscaled: DHFFC is flank-normalized, and the flanks carry the same
                    # ploidy as the variant.
                    if (vals[idx_ffc] + 0 < DEL_DHFFC + 0) supported = 1
                }
            }
            else if (idx_bfc > 0 && present(vals[idx_bfc])) {
                measured = 1
                # Scaled: DHBFC is normalized genome-wide, so a hemizygous baseline sits at
                # 0.5 and every threshold on it has to come down by the same factor.
                if (vals[idx_bfc] + 0 > DUP_DHBFC * copies / 2) supported = 1
            }
        }
        # Only tag where depth was actually measured for someone. No alignment means no
        # opinion, not a negative one.
        if (measured && !supported) {
            tags = (tags == "") ? "DEPTH_UNSUPPORTED" : tags ";DEPTH_UNSUPPORTED"
        }
    }

    $7 = (tags == "") ? "PASS" : tags
    print
}
