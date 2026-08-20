#!/usr/bin/env nextflow
nextflow.enable.types = true

/*
 * CADD + SpliceAI + SQUIRLS + ANNOVAR annotation of one cohort VCF per invocation.
 *
 * Input is a single already-merged VCF and its index, not a samplesheet -- this pipeline
 * annotates what it is given. To annotate more than one variant set, run it again with a
 * separate --vcf/--outdir per set.
 *
 * Rewritten from wf (~/git/wf/main.nf + modules/), whose annovar step ran all six ANNOVAR
 * commands as one monolithic per-interval process. Structure and conventions follow
 * the upstream pipeline this was ported from.
 */

include { mergeVersions } from '../shared/modules/versions.nf'
include { somalier_extract; somalier_extract_vcf; somalier_relate; somalier_ancestry;
          check_somalier } from '../shared/modules/qc/somalier.nf'
include { summary_report } from './modules/other/report.nf'
include { annotate_snps } from './subworkflows/annotate_snps.nf'

workflow {
    main:
    if (!params.vcf)       error "params.vcf (cohort VCF to annotate) is required"
    if (!params.tbi)       error "params.tbi (tabix index for params.vcf) is required"
    if (!params.cohort)    error "params.cohort is required -- it names the output files"
    if (!params.data_type) error "params.data_type is required -- it names the output files"

    // Site-local resources with no public default (conf/params.config leaves them null).
    // Checked up front, by name, so a missing one fails immediately instead of surfacing
    // later as a checkIfExists on a path the user never set -- but only when the stage that
    // reads them actually runs. Requiring a resource for work that will not happen is not
    // caution, it is a wrong error message.
    def required = [
        annovar_dir:                  "ANNOVAR install directory, containing humandb/ and the perl scripts",
        annovar_splice_scores_script: "annovar_parse_out_spliceai_squirls_scores.py (local path or pinned https:// URL)",
        omim_xref:                    "OMIM gene xref table for table_annovar.pl -xref",
    ]
    // A prescored table alone does not lift this: CADD still runs for whatever it does not
    // cover. Only the ANNOVAR-only re-run (table + skip flag) scores nothing.
    if (!(params.precomputed_cadd && params.skip_spliceai_squirls)) {
        required += [
            cadd_data_dir:            "CADD v1.6 GRCh38 data directory (must contain CADD.sh)",
        ]
    }
    if (!params.skip_spliceai_squirls) {
        required += [
            squirls_config:              "SQUIRLS config YAML",
            squirls_jannovar_model:      "SQUIRLS Jannovar transcript model (.ser)",
            spliceai_precomputed_scores: "SpliceAI precomputed scores VCF",
            spliceai_precomputed_tbi:    "tabix index for params.spliceai_precomputed_scores",
        ]
    }
    required.each { name, what ->
        if (!params[name]) error "params.${name} is required and has no default -- ${what}"
    }

    // Re-annotation entry points. Both halves of the CADD pair are needed or neither:
    // bcftools reads the table through its index, so a table without one fails inside a
    // task rather than at startup.
    if (params.precomputed_cadd as Boolean != params.precomputed_cadd_tbi as Boolean) {
        error "params.precomputed_cadd and params.precomputed_cadd_tbi must be given together"
    }

    // Sample QC (somalier), opt-in. The PED is OPTIONAL, and 04_qc gains ploidy.tsv and
    // cohort.inferred.ped either way -- what differs is the thing checked against: with a
    // PED, relatedness and sex are compared to what the pedigree claims; without one,
    // somalier's own inferred pedigree (--infer, via the draft ped) is cross-checked against
    // depth instead. Requiring a pedigree in order to look at the data would gate the check
    // on paperwork.
    snp_ancestry_params = ['somalier_labels', 'somalier_1kg_dir']
    snp_ancestry_configured = snp_ancestry_params.findAll { name -> params[name] }
    if (snp_ancestry_configured
            && snp_ancestry_configured.size() != snp_ancestry_params.size()) {
        error "somalier ancestry is partly configured: set both of " +
              "${snp_ancestry_params.join(', ')} or neither. Missing: " +
              "${(snp_ancestry_params - snp_ancestry_configured).join(', ')}"
    }
    if (snp_ancestry_configured && !params.somalier_sites) {
        error "params.somalier_sites is required for somalier ancestry -- ancestry is " +
              "predicted from the same extracted sites as relatedness"
    }
    // Alignments sidecar for somalier: when given, depth is read off the CRAM/BAMs instead
    // of FORMAT/AD in --vcf (see conf/params.config on why exomes need this). The reference
    // pair travels with the sheet, and a sheet without the sites file would silently do
    // nothing -- all three misuses are rejected at startup, by name.
    if (params.alignments && !params.alignment_reference) {
        error "params.alignment_reference is required when params.alignments is set -- " +
              "CRAM decode needs the exact FASTA the alignments were written against"
    }
    if (params.alignments && !params.alignment_reference_index) {
        error "params.alignment_reference_index is required when params.alignments is set"
    }
    if (params.alignments && !params.somalier_sites) {
        error "params.alignments does nothing without params.somalier_sites -- the somalier " +
              "QC stage is gated on the sites file"
    }
    // Nothing can verify the claim this flag makes -- that --vcf already carries SpliceAI and
    // SQUIRLS INFO. A plain VCF passed with it set yields empty splice columns, not an error,
    // so say plainly what was assumed rather than letting it pass in silence.
    if (params.skip_spliceai_squirls) {
        log.warn "skip_spliceai_squirls: assuming --vcf already carries SpliceAI and SQUIRLS " +
                 "INFO fields. Nothing checks this; if it does not, splice columns come out empty."
    }
    // The other unverifiable pairing: the flag together with a prescored table turns top-up
    // into passthrough. A table in CADD's own naming, or one missing variants, produces blank
    // CADD columns here rather than an error -- without the flag both are handled.
    if (params.precomputed_cadd && params.skip_spliceai_squirls) {
        log.warn "precomputed_cadd + skip_spliceai_squirls: the table is used as-is -- " +
                 "variants it does not cover get no CADD score, and its contigs must use " +
                 "the callset's naming. Nothing checks either. Without the flag, uncovered " +
                 "variants would be scored."
    }

    reference = record(
        fa:   file(params.reference.fa, checkIfExists: true),
        fai:  file(params.reference.fai, checkIfExists: true),
        dict: file(params.reference.dict, checkIfExists: true)
    )

    vcf                          = file(params.vcf, checkIfExists: true)
    tbi                          = file(params.tbi, checkIfExists: true)
    // null when the stage that would read it is skipped -- the guard block above has already
    // established that each is set whenever its stage runs. A staged path is only ever
    // dereferenced inside a process script, so a null that no process receives is never seen.
    cadd_data_dir                = params.cadd_data_dir
        ? file(params.cadd_data_dir, checkIfExists: true) : null
    squirls_config               = params.squirls_config
        ? file(params.squirls_config, checkIfExists: true) : null
    squirls_jannovar_model       = params.squirls_jannovar_model
        ? file(params.squirls_jannovar_model, checkIfExists: true) : null
    spliceai_precomputed_scores  = params.spliceai_precomputed_scores
        ? file(params.spliceai_precomputed_scores, checkIfExists: true) : null
    spliceai_precomputed_tbi     = params.spliceai_precomputed_tbi
        ? file(params.spliceai_precomputed_tbi, checkIfExists: true) : null
    annovar_dir                  = file(params.annovar_dir, checkIfExists: true)
    annovar_splice_scores_script = file(params.annovar_splice_scores_script, checkIfExists: true)
    omim_xref                    = file(params.omim_xref, checkIfExists: true)

    // --- ANNOVAR preflight -------------------------------------------------------------
    //
    // table_annovar.pl dies on a missing database AFTER CADD and SpliceAI have run -- hours
    // into a WGS run, for a fault visible in a directory listing at startup. Checked here,
    // in the workflow body: `file()` returns a java.nio.file.Path that resolves on the head
    // node, so the directory can be listed before any task is submitted. Not in a closure --
    // error() raised inside a .map{} surfaces as InvocationTargetException with the message
    // swallowed (measured; see annotate_svs/subworkflows/prepare_svs.nf).
    annovar_protocols = "${params.annovar_protocols}".replace('CLINVAR', "${params.clinvar_date}")

    def n_protocols = annovar_protocols.tokenize(',').size()
    def n_operations = "${params.annovar_operations}".tokenize(',').size()
    if (n_protocols != n_operations) {
        error "params.annovar_protocols has ${n_protocols} entries but params.annovar_operations " +
              "has ${n_operations} -- table_annovar.pl pairs them by position, so a mismatch " +
              "silently applies the wrong operation to every database after the first gap"
    }

    // ANNOVAR's filter databases fall back to a linear scan without a .idx, so .idx is NOT
    // required per protocol -- demanding it would be stricter than the tool and would fail a
    // working install. Two protocols break the hg38_<name>.txt convention and one is absent
    // by design; everything else follows it.
    def humandb = annovar_dir.resolve('humandb')
    def special = [
        // annotate_variation.pl maps the 1000g dbtypes onto this filename. A hg38_1000g*.txt
        // check fails on a correctly populated install.
        '1000g2015aug_all': ['hg38_ALL.sites.2015_08.txt'],
        // gx also needs the transcript FASTA.
        'refGene':          ['hg38_refGene.txt', 'hg38_refGeneMrna.fa'],
        // Built at runtime by build_cadd_humandb and symlinked in by stage_humandb.
        'CADDv1.6':         [],
    ]

    def missing = []
    annovar_protocols.tokenize(',').each { p ->
        (special.containsKey(p) ? special[p] : ["hg38_${p}.txt"]).each { f ->
            if (!humandb.resolve(f as String).exists()) missing << "${p} -> humandb/${f}"
        }
    }
    // stage_humandb symlinks *.txt, *.txt.idx and *.fa; an unmatched glob under
    // `set -euo pipefail` leaves the literal word and ln dies with a bare "No such file or
    // directory" deep in a task. At least one of each has to exist, even though no individual
    // .idx is required above.
    def humandb_names = humandb.list() as List<String>
    ['*.txt', '*.txt.idx', '*.fa'].each { g ->
        def re = g.replace('.', '\\.').replace('*', '.*')
        if (!humandb_names.any { n -> n ==~ re }) {
            missing << "no ${g} in humandb/ -- stage_humandb's symlink would fail"
        }
    }
    // Checked here rather than found late: only index_annovar.pl fails early, because
    // build_cadd_humandb is the first ANNOVAR-touching process.
    ['table_annovar.pl', 'variants_reduction.pl', 'convert2annovar.pl', 'index_annovar.pl'].each { s ->
        if (!annovar_dir.resolve(s).exists()) missing << "${s} not in ${params.annovar_dir}"
    }

    if (missing) {
        error "ANNOVAR preflight failed -- ${missing.size()} missing:\n  " +
              missing.join('\n  ') +
              "\n\nAll of them are listed at once so they can be fixed in one pass. " +
              "Note --clinvar_date is '${params.clinvar_date}', which names a database."
    }

    // [] when absent -- see the subworkflow's take: block on why a list, not a tuple.
    precomputed_cadd = params.precomputed_cadd
        ? [ file(params.precomputed_cadd, checkIfExists: true),
            file(params.precomputed_cadd_tbi, checkIfExists: true) ]
        : []

    // One item, so every process taking it as a singleton broadcasts against the
    // per-interval channels rather than consuming one interval's worth and stopping.
    vcf_tbi = channel.of(tuple(vcf, tbi))

    annotated = annotate_snps(
        vcf_tbi,
        reference, cadd_data_dir, squirls_config, squirls_jannovar_model,
        spliceai_precomputed_scores, spliceai_precomputed_tbi,
        "${params.cohort}", "${params.data_type}", "${params.clinvar_date}",
        annovar_protocols, "${params.annovar_operations}",
        omim_xref, annovar_dir, annovar_splice_scores_script,
        precomputed_cadd, params.skip_spliceai_squirls as Boolean
    )

    // Sample QC, off the cohort VCF -- or off the --alignments CRAM/BAMs when a sheet was
    // given, which is what makes the sex call usable on exomes. Independent of the
    // annotation stages -- nothing here feeds them, and they do not feed this -- so it is
    // gated on its own param and simply does not run without it.
    //
    // The processes are shared with annotate_svs (nf/shared/modules/qc); only the extract
    // step differs between the two pipelines. Called inline rather than through a
    // subworkflow: it is a few calls and two conditionals.
    if (params.somalier_sites) {
        // [] when no PED was given -- see somalier_relate's `ped` input on why a list.
        qc_ped = params.ped ? [ file(params.ped, checkIfExists: true) ] : []
        qc_sites = file(params.somalier_sites, checkIfExists: true)

        // Only the extract step branches; everything below consumes .somalier files and
        // does not care where they came from.
        if (params.alignments) {
            // One task per sheet row. The sheet columns are the contract shared with
            // annotate_svs' --alignments (sample,alignment,alignment_index); parsed eagerly,
            // not in a .map{}, because error() inside an operator closure surfaces as
            // InvocationTargetException with the message swallowed.
            def aln_rows = file(params.alignments, checkIfExists: true).splitCsv(header: true)
            if (!aln_rows) error "params.alignments has no data rows: ${params.alignments}"
            aln_rows.eachWithIndex { row, i ->
                ['sample', 'alignment', 'alignment_index'].each { col ->
                    if (!row[col]) {
                        error "--alignments row ${i + 1} is missing '${col}' -- the sheet " +
                              "needs sample,alignment,alignment_index"
                    }
                }
            }
            aln_parts = channel.fromList(aln_rows.collect { row ->
                tuple((row.sample as String).trim(),
                      file((row.alignment as String).trim(), checkIfExists: true),
                      file((row.alignment_index as String).trim(), checkIfExists: true))
            }).multiMap { s, a, x ->
                sample: s
                alignment: a
                index: x
            }
            extracted = somalier_extract(
                aln_parts.sample, aln_parts.alignment, aln_parts.index, qc_sites,
                file(params.alignment_reference, checkIfExists: true),
                file(params.alignment_reference_index, checkIfExists: true))
            // Sorted before it reaches a task: Nextflow emits collected items in completion
            // order, and an order that varies between runs defeats -resume for no benefit.
            qc_files = extracted.out.toSortedList { a, b -> a.name <=> b.name }
            // Every per-sample extract task reports the same build; one copy is enough.
            qc_extract_versions = extracted.versions.first()
        }
        else {
            extracted = somalier_extract_vcf(vcf, tbi, qc_sites, reference.fa, reference.fai)
            // Sorted for the same -resume reason; this extract already emits one List<Path>.
            qc_files = extracted.out.map { fs -> (fs as List).sort { a, b -> a.name <=> b.name } }
            qc_extract_versions = extracted.versions
        }
        related = somalier_relate("${params.cohort}", qc_files, qc_ped, !params.ped)

        // check_somalier always runs: with an operator PED it is the pedigree COMPARISON --
        // ploidy against the PED's sex column, and an inferred PED to replace a contradicted
        // one. Without one, somalier's own inferred pedigree (--infer, via the draft ped) is
        // cross-checked against depth instead, so 04_qc gains a ploidy.tsv and a
        // cohort.inferred.ped the SV pipeline's --ped can consume either way.
        check_somalier_awk = file("${projectDir}/../shared/assets/check_somalier.awk",
                                  checkIfExists: true)
        checked = check_somalier(related.samples, related.pairs,
                                 qc_ped ? qc_ped[0] : related.draft_ped, check_somalier_awk)
        qc_ploidy_out = checked.ploidy
        qc_inferred_ped_out = checked.inferred_ped
        qc_check_versions = checked.versions

        if (snp_ancestry_configured) {
            ancestry = somalier_ancestry(
                "${params.cohort}",
                qc_files,
                file(params.somalier_labels, checkIfExists: true),
                file(params.somalier_1kg_dir, checkIfExists: true))
            qc_ancestry_out = ancestry.out
            qc_ancestry_cohort = ancestry.cohort_out
            qc_ancestry_page = ancestry.html
            qc_ancestry_versions = ancestry.versions
        }
        else {
            qc_ancestry_out = channel.empty()
            qc_ancestry_cohort = channel.empty()
            qc_ancestry_page = channel.empty()
            qc_ancestry_versions = channel.empty()
        }

        qc_samples_out = related.samples
        qc_pairs_out = related.pairs
        qc_html_out = related.html
        qc_versions = qc_extract_versions
            .mix(related.versions, qc_check_versions, qc_ancestry_versions)
    }
    else {
        qc_ploidy_out = channel.empty()
        qc_inferred_ped_out = channel.empty()
        qc_samples_out = channel.empty()
        qc_pairs_out = channel.empty()
        qc_html_out = channel.empty()
        qc_ancestry_out = channel.empty()
        qc_ancestry_cohort = channel.empty()
        qc_ancestry_page = channel.empty()
        qc_versions = channel.empty()
    }

    // The run summary: one tidy TSV over the input VCF, the multianno TSV and whatever QC
    // tables exist. The somalier channels feed it via toList(), which turns "the QC stage
    // did not run" (an empty channel) into the empty list its optional inputs expect --
    // no re-branching on params here.
    def run_date = "${params.clinvar_date}".tokenize('_')[1]
    report = summary_report(
        vcf, tbi, annotated.annovar_tsv,
        qc_samples_out.toList(), qc_ploidy_out.toList(),
        qc_ancestry_cohort.toList(), qc_pairs_out.toList(),
        "${params.cohort}", "${params.data_type}", run_date, "${params.clinvar_date}")

    // One versions.yml per distinct process, collected from every stage and published to
    // pipeline_info/. mergeVersions deduplicates, so a process that ran 24 per-interval
    // tasks appears once.
    software_versions = mergeVersions(annotated.versions.mix(qc_versions, report.versions))

    publish:
    // cadd_merged is one (tsv.gz, tbi) tuple; split so each lands as its own published
    // file rather than a nested collection.
    cadd_tsv          = annotated.cadd_merged.map { t, _i -> t }
    cadd_tbi          = annotated.cadd_merged.map { _t, i -> i }
    annotations       = annotated.annovar_tsv
    // Numbered 03 despite being produced before the ANNOVAR stage: renumbering 02_annovar
    // would break the path the README and the plan's acceptance criteria both name.
    annotated_vcf     = annotated.annotated_vcf.map { v, _i -> v }
    annotated_vcf_tbi = annotated.annotated_vcf.map { _v, i -> i }
    qc_ploidy         = qc_ploidy_out
    qc_inferred_ped   = qc_inferred_ped_out
    qc_samples_tsv    = qc_samples_out
    qc_pairs_tsv      = qc_pairs_out
    qc_relate_html    = qc_html_out
    // The full table (this run's samples plus the 1kg reference set) and the same table cut
    // down to this run. Both are published: the cohort one is what you read, the full one is
    // the background that makes it interpretable.
    qc_ancestry_tsv   = qc_ancestry_out
    qc_ancestry_cohort_tsv = qc_ancestry_cohort
    // Not named qc_ancestry_page: a publish target sharing a name with the local it reads
    // resolves to itself.
    qc_ancestry_html  = qc_ancestry_page
    summary_tsv       = report.out
    pipeline_versions = software_versions
}
output {
    cadd_tsv          { path "01_cadd" }
    cadd_tbi          { path "01_cadd" }
    annotations       { path "02_annovar" }
    annotated_vcf     { path "03_annotated_vcf" }
    annotated_vcf_tbi { path "03_annotated_vcf" }
    // Sample QC, present only when --somalier_sites was given. Nothing downstream in this
    // pipeline reads it; it is published because a sample swap or a wrong pedigree is
    // invisible in an annotated VCF.
    qc_ploidy         { path "04_qc" }
    qc_inferred_ped   { path "04_qc" }
    qc_samples_tsv    { path "04_qc" }
    qc_pairs_tsv      { path "04_qc" }
    qc_relate_html    { path "04_qc" }
    qc_ancestry_tsv   { path "04_qc" }
    qc_ancestry_cohort_tsv { path "04_qc" }
    qc_ancestry_html  { path "04_qc" }
    summary_tsv       { path "05_report" }
    pipeline_versions { path "pipeline_info" }
}
