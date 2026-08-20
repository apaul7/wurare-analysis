nextflow.enable.types = true
include { GroupedVcf } from '../modules/prepare/standardize.nf'
// One process definition, two call sites. Nextflow refuses to invoke a process twice in
// one workflow ("has been already used"), so the two axes take aliases -- the alias is the
// requirement, not a second copy of the logic.
include { MergedVcf; fill_tags; passthrough } from '../modules/merge/svdb_merge.nf'
include { svdb_merge as merge_axis_a; svdb_merge as merge_axis_b } from '../modules/merge/svdb_merge.nf'
include { promote_support as promote_axis_a; promote_support as promote_axis_b } from '../modules/merge/svdb_merge.nf'

// Phase 2: axis A then axis B.
//
// Two axes and not one. Merging both in a single pass makes a two-caller/one-sample
// duplicate indistinguishable from a one-caller/two-sample shared event -- axis A asks "is
// this call real", axis B asks "which samples carry it".
//
// The two stages are deliberately decoupled: axis B is one process invocation over one
// channel, so swapping it for `gatk SVCluster` later touches this file and nothing else.
workflow merge_svs {
    take:
    standardized: Channel<GroupedVcf>
    promote_awk: Path
    overlap_axis_a: String
    overlap_axis_b: String
    bnd_distance: String
    ins_distance: String

    main:
    // Group on the SET of sample IDs, not sample_set and not sample count. Two
    // callers meet at axis A when their sample sets are identical; anything else goes
    // straight to axis B.
    //
    // Sorted by caller inside each group, and the groups themselves sorted, because
    // Nextflow emits collected items in completion order -- unsorted, two runs over
    // identical inputs can produce different cohort VCFs.
    grouped = standardized
        .map { g -> tuple(g.sample_key, g.entry) }
        .groupTuple()
        .map { key, entries ->
            // Axis-A --priority: the local trust order. This is
            // not a tie-break. SVDB gives the merged record the priority input's
            // representation -- its breakpoints and, for a sample both inputs genotyped, its
            // genotype (measured during design) -- so the order is a claim about whose call to
            // believe, and it is site-local knowledge rather than a property of the tools.
            // Hardcoded here for now; making it a param is an open item.
            //
            // A caller absent from this list sorts after every caller in it, alphabetically
            // among its peers, so adding one to the sheet can never silently outrank a
            // trusted caller -- it can only land behind them.
            // Matched case- and whitespace-insensitively: a casing difference is a typo, not
            // a different tool. An exact match sorted "Manta" behind every other caller.
            def trusted = ['manta', 'smoove', 'cnvkit', 'cnvnator']
            def ordered = (entries as List).sort { a, b ->
                def ra = trusted.indexOf((a.caller as String).trim().toLowerCase())
                def rb = trusted.indexOf((b.caller as String).trim().toLowerCase())
                ra = ra < 0 ? trusted.size() : ra
                rb = rb < 0 ? trusted.size() : rb
                ra == rb ? (a.caller as String) <=> (b.caller as String) : ra <=> rb
            }
            // The label names files and becomes the axis-B tag, so it MUST be unique per
            // group. sample_set alone is not: two different sample sets can carry the same
            // free-text label, and then both groups emit a file called <label>.vcf.gz, one
            // silently overwrites the other in the axis-B work dir, and SVDB is handed the
            // same file twice under a duplicated tag. Suffixing with a short digest of the
            // sample_key keeps the label readable and makes the collision impossible.
            //
            // It also has to be a legal VCF tag name, because axis B hands it to SVDB as the
            // `--vcf file:tag` tag and SVDB prefixes real INFO keys with it (<tag>_CHROM,
            // <tag>_POS, <tag>_INFO ...). htslib takes only letters, digits and underscore
            // there and will not take a leading digit, so a sample_set like "sample.2"
            // yields "sample.2_<hash>_INFO" and bcftools rejects the header line when
            // promote_support reads the cohort VCF back. The digest is over the sample_key,
            // not over this stem, so folding characters together cannot merge two groups.
            def suffix = (key as String).md5().substring(0, 6)
            // Named from the alphabetically first sample_set, NOT from ordered[0] -- the
            // label ends up in filenames and in the axis-B tag, so tying it to the priority
            // order would rename every output the moment a caller is added to the sheet.
            def named = (entries as List).collect { e -> e.sample_set as String }.min()
            def raw = named.replaceAll(/[^A-Za-z0-9_]/, '_')
            def stem = raw ==~ /^[0-9].*/ ? "_${raw}" : raw
            // A group is joint if any input in it came from a joint caller; that decides
            // where it sorts at axis B, below.
            def joint = (entries as List).any { e -> e.joint }
            tuple("${stem}_${suffix}" as String, ordered, joint)
        }

    forked = grouped.branch { _key, entries, _joint ->
        merge: (entries as List).size() > 1
        single: true
    }

    // Axis A: cross-caller reconciliation. Tag is the caller; priority is the sort order
    // established above -- manta, then the rest.
    axis_a_in = forked.merge.map { key, entries, joint ->
        tuple(key,
              (entries as List).collect { e -> e.vcf },
              (entries as List).collect { e -> e.caller as String },
              joint)
    }
    axis_a = merge_axis_a(
        axis_a_in.map { key, _v, _t, _j -> key },
        axis_a_in.map { _k, vcfs, _t, _j -> vcfs },
        axis_a_in.map { _k, _v, tags, _j -> tags },
        overlap_axis_a, bnd_distance, ins_distance
    )

    // Mandatory between the axes (verified by experiment during design): axis B appends
    // its own svdb_origin/FOUNDBY/SUPP_VEC/set/VARID without replacing axis A's, and
    // bcftools then returns the stale axis-A value. Promote to CALLER_SUPP/NCALLER and strip the raw keys first.
    // JOINED on the label, not zipped on arrival. axis_a.out emits in task COMPLETION order
    // while axis_a_in emits in channel order, so feeding both to one process positionally
    // pairs the nth label with the nth-FINISHED VCF. That made this the one step -resume
    // could never match -- its hash changed every run -- and it silently risked handing a
    // group's VCF another group's label and `joint`, which decides axis-B priority.
    // svdb_merge now carries its label on the item so the pairing can be keyed.
    axis_a_meta = axis_a_in.map { key, _v, _t, joint -> tuple(key, joint) }
    axis_a_paired = axis_a.out
        .map { m -> tuple(m.label, m.vcf) }
        .join(axis_a_meta, by: 0)   // typed syntax requires `by` to be explicit

    axis_a_parts = axis_a_paired.multiMap { key, vcf, joint ->
        label: key
        vcf: vcf
        joint: joint
    }
    promoted = promote_axis_a(
        axis_a_parts.label,
        axis_a_parts.vcf,
        promote_awk,
        channel.value("axis_a"),
        axis_a_parts.joint
    )

    // A lone input has nothing to reconcile; it is repacked so axis B's inputs are uniform.
    singles = passthrough(
        forked.single.map { key, _entries, _joint -> key },
        forked.single.map { _key, entries, _joint -> (entries as List)[0].vcf },
        forked.single.map { _key, _entries, joint -> joint }
    )

    // Axis B: cohort assembly, joint inputs first. This is the joint-caller priority
    // decision, and until now the code did not implement it: sorting on the label alone
    // handed priority to whichever sample_set sorted first, so a sample present in both a
    // joint VCF and its own single-sample calls kept the joint genotype only by luck of the
    // alphabet. The design spike measured that priority order changes the surviving genotype
    // (0/1 vs 1/1 on one fixture), which is the whole reason joint VCFs are kept unsplit.
    // Label breaks the tie so that two runs still agree.
    axis_b_inputs = promoted.out.mix(singles.out)
        .toSortedList { a, b ->
            a.joint == b.joint ? a.label <=> b.label : (a.joint ? -1 : 1)
        }

    cohort_merged = merge_axis_b(
        channel.value("cohort"),
        axis_b_inputs.map { entries -> (entries as List).collect { e -> e.vcf } },
        axis_b_inputs.map { entries -> (entries as List).collect { e -> e.label as String } },
        overlap_axis_b, bnd_distance, ins_distance
    )

    // MODE=cohort: axis B's svdb_origin names sample sets, so promoting it would overwrite
    // caller support with sample-set labels. The caller union comes from the per-input
    // blobs instead. See the header of promote_caller_support.awk.
    // joint=true past this point: the cohort VCF carries one genotype matrix over every
    // sample, whatever it was assembled from. Nothing downstream reads the flag -- it exists
    // for the axis-B sort above -- but false here would be the wrong answer to a question
    // someone may yet ask.
    // Safe to read positionally: axis B is ONE invocation over one collected channel, so
    // there is no arrival order to get wrong. Only the .vcf is needed -- the label here is
    // the constant "cohort", not the record's.
    cohort_promoted = promote_axis_b(channel.value("cohort"),
                                     cohort_merged.out.map { m -> m.vcf },
                                     promote_awk, channel.value("cohort"),
                                     channel.value(true))

    // Internal AC/AN/AF in the main line, not the Talos tail.
    //
    // Destructured rather than passed as a MergedVcf: a record reaching a process `input:`
    // is hashed by object identity, so its task hash changes every run and -resume can
    // never match it. Records stay in the channel plumbing; processes take primitives.
    // .first() because the four reads below would otherwise each try to consume the same
    // queue channel.
    promoted_cohort = cohort_promoted.out.first()
    cohort = fill_tags(
        promoted_cohort.map { c -> c.label },
        promoted_cohort.map { c -> c.vcf },
        promoted_cohort.map { c -> c.tbi },
        promoted_cohort.map { c -> c.joint }
    )

    all_versions = axis_a.versions
        .mix(promoted.versions, singles.versions, cohort_merged.versions,
             cohort_promoted.versions, cohort.versions)

    emit:
    cohort: Channel<MergedVcf> = cohort.out
    versions: Channel<Path> = all_versions
}
