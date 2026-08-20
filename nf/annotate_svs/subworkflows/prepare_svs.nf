nextflow.enable.types = true
include { SvInput; AlignmentEntry } from '../modules/other/samplesheet.nf'
include { StandardVcf; GroupedVcf; read_sample_ids; pass_size_filter; finalize } from '../modules/prepare/standardize.nf'

// Phase 1: every input VCF onto one schema, plus the preflight facts later phases are
// entitled to assume.
//
// No splitting of multi-sample VCFs: joint calls carry the inheritance signal, and
// splitting destroys the 0/0-vs-./. distinction that inheritance reasoning runs on.
workflow prepare_svs {
    take:
    sv_inputs: Channel<SvInput>
    alignments: Channel<AlignmentEntry>
    ped: Path
    normalize_awk: Path
    stamp_awk: Path
    min_sv_size: Integer

    main:
    // Sample IDs are the join key across every input. Read them from the headers
    // first, because the grouping the merge will use is derived from them, not from the
    // sheet's free-text sample_set label.
    // Destructured with multiMap, twice over. A record in a process `input:` block is hashed
    // by object identity, so its task hash changes every run and -resume never matches it --
    // records stay in the channel plumbing and in `output:`, never in an `input:`.
    //
    // multiMap rather than five separate .map calls: one traversal, and the five fields stay
    // visibly one row. (Five .map reads would also work -- a queue channel is broadcast to
    // every consumer, which is why `sv_inputs` can feed two processes below.) Forked twice
    // because a fork set is consumed by one process.
    rows_ids = sv_inputs.multiMap { i ->
        sample_set: i.sample_set
        caller: i.caller
        joint: i.joint
        vcf: i.vcf
        tbi: i.tbi
    }
    with_samples = read_sample_ids(rows_ids.sample_set, rows_ids.caller, rows_ids.joint,
                                   rows_ids.vcf, rows_ids.tbi)

    // Grouping key: the SET of sample IDs a VCF contains -- not sample count, not
    // sample_set. Two callers' VCFs meet at axis A when their sample sets are identical; a
    // VCF whose sample set matches nothing goes straight to axis B.
    keyed = with_samples.out.map { e ->
        def names = e.samples.text.readLines().collect { l -> l.trim() }.findAll { l -> l }
        tuple(names.sort().join('|'), names, e)
    }

    // Preflight is recorded, not merely computed: every failure mode here is "it happened
    // quietly", so the routing decisions are written down where a reader sees them.
    //
    // Collected into a published file rather than .view()'d onto stdout. At cohort scale the
    // grouping line names every sample in every group, and the ped and alignment lines name
    // every sample missing a row -- which buries the Nextflow progress output it is printed
    // among. A file is where a record of what the run decided belongs anyway: it is wanted
    // after the fact, not during.
    per_input = keyed.map { _key, names, input ->
        "preflight   ${input.sample_set}\t${input.caller}\tjoint=${input.joint}" +
        "\tsamples=${names.size()}" +
        "\n"
    }

    grouping = keyed.map { key, names, _input -> tuple(key, names) }
        .groupTuple()
        .map { _key, names ->
            "grouping    ${names.size()} input(s) share [${(names[0] as List).join(', ')}]\n"
        }

    // Samples with no alignment skip the depth stages. Listed, never dropped quietly.
    all_samples = keyed.map { _key, names, _input -> names }.collect(flat: true)
    aligned = alignments.map { a -> a.sample }.collect(flat: true).ifEmpty([])
    // unique(false), never unique(). Groovy's no-arg unique() mutates the list in place and
    // returns it, and `as List` on something already a List hands back the SAME object -- so
    // this and the ped check below were both mutating the one list that the `all_samples`
    // value channel gives to every consumer. Two closures mutating it at once is a
    // ConcurrentModificationException, and it only shows up once the cohort is big enough for
    // the two to overlap in time. Small runs got away with it.
    alignment_line = all_samples.combine(aligned).map { all, have ->
        def missing = (all as List).unique(false).findAll { s -> !(have as List).contains(s) }
        missing
            ? "alignments  missing for ${missing.sort().join(', ')} -- depth stages skip these\n"
            : "alignments  present for every sample\n"
    }

    // The other direction. A row naming a sample this cohort lacks is not an error -- one
    // site-wide alignments.csv reused across cohorts is normal -- but it must be dropped
    // before `bcftools view -s`, which hard-errors on an unknown sample. Named out loud so a
    // typo'd sample ID does not read as a sample that simply has no depth.
    //
    // TOTAL mismatch is the exception to that, and it is fatal. main.nf intersects the sheet
    // against this same list, so with no overlap `in_cohort` is empty: the depth stages then
    // run over nothing, and somalier_relate/_ancestry are handed an empty file list, which
    // somalier 0.2.19 dies on with "index out of bounds, the container is empty" -- its
    // ancestry.nim indexes s[0] to read a column count off a zero-length matrix, so the
    // traceback names neither the sheet nor the cohort. Checked here because this closure
    // already holds both lists.
    //
    // In a closure, not the workflow body: all_samples comes from a process, so unlike
    // ped_samples below it cannot be read eagerly. error() from HERE does surface its message
    // -- measured on 26.04.6, exit 1 with the text intact, with the closure firing after an
    // upstream process rather than at DAG construction. That is not in tension with the ped
    // block's note; it is a different construct, so do not "fix" this to match it.
    extra_line = all_samples.combine(aligned).map { all, have ->
        def cohort = (all as List).unique(false)
        def sheet = (have as List).unique(false)
        if (sheet && !sheet.any { s -> cohort.contains(s) }) {
            error "params.alignments: no sample in the sheet appears in this cohort, so " +
                  "every depth and somalier stage would run over nothing. " +
                  "Sheet: ${sheet.sort().join(', ')}. " +
                  "Cohort: ${cohort.sort().join(', ')}"
        }
        def extra = sheet.findAll { s -> !cohort.contains(s) }
        extra
            ? "alignments  ${extra.sort().join(', ')} are in the sheet but not in the cohort " +
              "-- ignored\n"
            : ""
    }

    // The PED is required, and a required file nothing reads is theatre. A sample with
    // no row has no sex, so it is invisible to every ploidy-aware and sex-stratified step
    // downstream -- silently, because those steps see a missing key rather than an error.
    // Reported here alongside the other routing facts. Parsed in the workflow body rather
    // than a process: it is a handful of lines, and error() raised inside a .map{} closure
    // surfaces as an InvocationTargetException with the real message swallowed.
    ped_samples = ped.text.readLines()
        .findAll { l -> l.trim() && !l.trim().startsWith('#') }
        .collect { l ->
            def fields = l.trim().split(/\s+/)
            if (fields.size() < 6) {
                error "params.ped: expected 6 whitespace-separated columns (family, " +
                      "sample, father, mother, sex, phenotype), got ${fields.size()} " +
                      "in line: ${l.trim()}"
            }
            fields[1]
        }

    ped_line = all_samples.map { all ->
        def missing = (all as List).unique(false).findAll { s -> !ped_samples.contains(s) }
        missing
            ? "ped         no row for ${missing.sort().join(', ')} -- no sex, so no " +
              "sex-stratified AF and no ploidy-aware handling for them\n"
            : "ped         ${ped_samples.size()} row(s), covering every cohort sample\n"
    }

    // sort: true so the file is deterministic and the four kinds of line group together,
    // rather than interleaving in whatever order the channels happened to emit.
    preflight_report = per_input.mix(grouping, alignment_line, extra_line, ped_line)
        .collectFile(name: 'preflight.txt', sort: true)

    // PASS,. + size floor before normalization, so the counts either side of the filter
    // are attributable to one input.
    rows_filter = sv_inputs.multiMap { i ->
        sample_set: i.sample_set
        caller: i.caller
        joint: i.joint
        vcf: i.vcf
        tbi: i.tbi
    }
    filtered = pass_size_filter(rows_filter.sample_set, rows_filter.caller,
                                rows_filter.joint, rows_filter.vcf, rows_filter.tbi,
                                min_sv_size)

    // One path for every caller. There is no svtk branch any more -- svtk standardize is
    // unrunnable, so normalize_records.awk is the only standardization and it applies
    // uniformly. Fewer paths, and no branch that silently never fires.
    fin = filtered.out.multiMap { e ->
        sample_set: e.sample_set
        caller: e.caller
        joint: e.joint
        filtered_vcf: e.vcf
    }
    ready = finalize(fin.sample_set, fin.caller, fin.joint, fin.filtered_vcf,
                     normalize_awk, stamp_awk)

    // Collected in main: rather than inline in emit: -- a chained .mix() in an emit
    // expression does not resolve under typed syntax ("No such variable").
    all_versions = with_samples.versions
        .mix(filtered.versions, ready.versions)

    // Attach the sample-ID key to each standardized VCF. Until this existed the key was
    // computed here, printed, and then thrown away -- while merge_svs grouped on the sheet's
    // sample_set label instead. Two callers on the same samples labelled differently then
    // never met at axis A, and worse, two DIFFERENT sample sets sharing a label were merged
    // across at cross-caller thresholds. Joined on (sample_set, caller), which identifies a
    // sheet row uniquely.
    key_by_row = keyed.map { key, _names, input ->
        tuple("${input.sample_set}\u0000${input.caller}" as String, key)
    }
    grouped_out = ready.out
        .map { e -> tuple("${e.sample_set}\u0000${e.caller}" as String, e) }
        .join(key_by_row, by: 0)   // typed syntax requires `by` to be explicit
        .map { _row, e, key -> new GroupedVcf(sample_key: key as String, entry: e) }

    emit:
    standardized: Channel<GroupedVcf> = grouped_out
    // Emitted so main can intersect the alignments sheet against it -- see extra_line above.
    cohort_samples: Channel<List<String>> = all_samples
    filter_counts: Channel<Path> = filtered.counts
    preflight: Channel<Path> = preflight_report
    versions: Channel<Path> = all_versions
}
