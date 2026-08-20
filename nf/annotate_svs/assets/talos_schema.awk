# Rename and derive the INFO fields Talos requires.
#
# Usage:  awk -v GNOMAD_POP=gnomad_v4.1 -f talos_schema.awk cohort.vcf
#
# Talos reads a fixed set of field names (traced from run_hail_filtering_sv.py). This pipeline
# already computes the same quantities under its own names, so the tail is a rename plus two
# derivations -- not a re-annotation:
#
#   gnomad_sv_AF        ->  {GNOMAD_POP}_sv_AF     population frequency, hard filter < 0.03
#   (from genotypes)    ->  N_HET, N_HOMALT        carrier counts
#   (stamped if absent) ->  STATUS
#   FILTER tags         ->  SOFT_FILTERS, FILTER=PASS   Talos deletes a non-empty FILTER
#
# AC/AN/AF already come from bcftools +fill-tags at Phase 2, and ALGORITHMS from Phase 1 --
# which matters more than it looks, because Talos DEFAULTS ALGORITHMS to ['gCNV'] when absent.
# An unstamped callset is therefore mislabelled rather than rejected.
#
# N_HET and N_HOMALT are counted here rather than taken from `bcftools +fill-tags -t AC_Het`
# because those count ALLELES, not samples: AC_Hom is twice the homozygote count. Talos wants
# sample counts, and a factor-of-two error in a carrier count is exactly the kind of thing
# that reads as plausible and is wrong.
#
# NOT done here: PREDICTED_LOF, which only gatk SVAnnotate can produce and which Talos hard
# gates on -- a variant with no PREDICTED_LOF entry is dropped entirely. That runs upstream of
# this in the same subworkflow. Nor the AF_MALE/AF_FEMALE *values*, which need sample sex from
# a PED and come from bcftools +fill-tags upstream -- only their headers are backfilled here.
#
# Header declaration is the real contract, read off Talos's rearrange_annotations(): it
# tolerates a missing ALGORITHMS, STATUS, CHR2 and END2, and nothing else. Every other field
# is read as `mt.info.X` / `mt.info[f'{GNOMAD_POP}_sv_SVID']`, and hail builds that struct from
# the VCF *header*, not from the records. An undeclared field is a hail error at annotation
# time -- Talos does not start. Absent per-record VALUES are fine, hail reads them as missing.
# So declaring a field this pipeline cannot populate is strictly better than omitting it: it
# turns a crash into a null.
#
# Two fields are declared without values for exactly that reason:
#   {GNOMAD_POP}_sv_SVID  svdb query returns an occurrence count and a frequency, not the
#                         identifier of the matching gnomAD variant. Talos reports SVID rather
#                         than filtering on it, so a null is a degraded report -- an
#                         undeclared one is no report at all.
#   AF_MALE / AF_FEMALE   fill-tags declares only the groups the PED actually contains. A
#                         single-sex cohort leaves the other undeclared, and Talos branches on
#                         seeing AF_MALE, then reads AF_FEMALE off the same struct. One sex
#                         present is the crash case, not zero.

# Declare an INFO field only if the input has not already declared it. Emitting it twice is
# its own bug: bcftools keeps the first and warns, so a second declaration with a different
# Number/Type is silently ignored rather than rejected.
function declare(id, number, type, description) {
    if (id in declared) return
    print "##INFO=<ID=" id ",Number=" number ",Type=" type ",Description=\"" description "\">"
    declared[id] = 1
}

function info_has(info, key,    n, i, parts) {
    n = split(info, parts, ";")
    for (i = 1; i <= n; i++) if (parts[i] ~ ("^" key "=")) return 1
    return 0
}

function info_get(info, key,    n, i, parts) {
    n = split(info, parts, ";")
    for (i = 1; i <= n; i++) if (parts[i] ~ ("^" key "=")) return substr(parts[i], length(key) + 2)
    return ""
}

BEGIN {
    FS = OFS = "\t"
    if (GNOMAD_POP == "") GNOMAD_POP = "gnomad_v4.1"
}

/^##INFO=<ID=/ {
    id = $0
    sub(/^##INFO=<ID=/, "", id)
    sub(/[,>].*$/, "", id)
    declared[id] = 1
    print
    next
}

/^#CHROM/ {
    declare(GNOMAD_POP "_sv_AF", 1, "Float", "Population AF, Talos schema")
    declare(GNOMAD_POP "_sv_SVID", 1, "String", \
            "Matching gnomAD-SV variant ID, Talos schema -- declared unpopulated, see talos_schema.awk")
    # Declared defensively, not because this stage writes them. normalize_records.awk derives
    # SVTYPE/SVLEN/END onto records whose caller may never have declared them, and a caller
    # header is the one input here that no process in this repo produces. AC/AN/AF are left
    # out on purpose: fill_tags writes those unconditionally at Phase 2, so backfilling them
    # would only hide it having not run.
    declare("SVTYPE", 1, "String", "Type of structural variant")
    declare("SVLEN", 1, "Integer", "Length of structural variant")
    declare("END", 1, "Integer", "End position of structural variant")
    declare("SOFT_FILTERS", ".", "String",
            "FILTER tags moved here for Talos, which drops any record with a non-empty FILTER")
    declare("N_HET", 1, "Integer", "Number of heterozygous samples")
    declare("N_HOMALT", 1, "Integer", "Number of homozygous-alt samples")
    declare("STATUS", 1, "String", "Call status")
    # Number=A to match what bcftools +fill-tags declares when it does compute these, so the
    # backfilled header and the computed one never disagree on arity.
    declare("AF_MALE", "A", "Float", "Allele frequency in male samples")
    declare("AF_FEMALE", "A", "Float", "Allele frequency in female samples")
    print
    next
}

/^#/ { print; next }

{
    info = ($8 == ".") ? "" : $8

    # FILTER to INFO, then PASS. Talos runs
    #     mt.filter_rows(hl.is_missing(mt.filters) | (mt.filters.length() == 0))
    # BEFORE any category logic, so a soft tag is not soft there -- it is deletion. The pipeline puts
    # tags in FILTER precisely because nothing should ever be removed, and that promise
    # inverts at this one consumer: a fully annotated callset reaches Talos and vanishes,
    # with nothing anywhere reporting it. Measured, not theorised -- a 2-record VCF tagged
    # COMMON_INTERNAL loaded as 0 rows.
    #
    # The tags are carried, not discarded: FILTER is where every other consumer reads them,
    # and this is the Talos-specific translation stage, so the move happens here and nowhere
    # upstream. Semicolons become commas because a semicolon separates INFO fields.
    if ($7 != "" && $7 != "." && $7 != "PASS") {
        if (!info_has(info, "SOFT_FILTERS")) {
            soft = $7
            gsub(/;/, ",", soft)
            info = info ";SOFT_FILTERS=" soft
        }
        $7 = "PASS"
    }

    # Population AF under the name Talos looks for. Copied rather than moved: the pipeline's
    # own consumers still read gnomad_sv_AF, and two names for one number is cheaper than a
    # rename that breaks the other half of the output.
    pop = info_get(info, "gnomad_sv_AF")
    if (pop != "" && pop != "." && !info_has(info, GNOMAD_POP "_sv_AF")) {
        info = info ";" GNOMAD_POP "_sv_AF=" pop
    }

    # Carrier counts, from the genotypes themselves.
    n = split($9, keys, ":")
    gt_idx = 0
    for (i = 1; i <= n; i++) if (keys[i] == "GT") gt_idx = i

    n_het = 0; n_homalt = 0
    if (gt_idx > 0) {
        for (c = 10; c <= NF; c++) {
            split($c, vals, ":")
            gt = vals[gt_idx]
            gsub(/\|/, "/", gt)
            if (gt == "0/1" || gt == "1/0") n_het += 1
            else if (gt == "1/1") n_homalt += 1
        }
    }
    if (!info_has(info, "N_HET"))    info = info ";N_HET=" n_het
    if (!info_has(info, "N_HOMALT")) info = info ";N_HOMALT=" n_homalt

    if (!info_has(info, "STATUS")) info = info ";STATUS=PASS"

    sub(/^;/, "", info)
    $8 = (info == "") ? "." : info
    print
}
