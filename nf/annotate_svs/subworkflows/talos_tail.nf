nextflow.enable.types = true
include { MergedVcf } from '../modules/merge/svdb_merge.nf'
include { svannotate; sex_stratified_af; talos_schema } from '../modules/talos/talos_tail.nf'

// Phase 6: Talos-schema output, optional and separate.
//
// Talos is ONE consumer, not the goal. Everything upstream serves every consumer; this
// converts the shared cohort VCF into the only schema Talos accepts. Skipping it costs
// nothing else.
workflow talos_tail {
    take:
    cohort: Channel<MergedVcf>
    protein_coding_gtf: Path
    noncoding_bed: Path
    // A channel, not a Path: main.nf passes the PED carrying somalier's inferred sex, which
    // is produced by a process rather than staged from a param. Identical in content to
    // params.ped when somalier did not run.
    ped: Channel<Path>
    symbolic_awk: Path
    schema_awk: Path
    keep_awk: Path
    check_awk: Path
    gnomad_pop: String

    main:
    // Value channel: `cohort` is read three times, and a queue channel can only be consumed
    // once. Same defect as annotate_cohort had, and for the same reason -- this path needs a
    // GTF and noncoding BED, so it had never run.
    coh = cohort.first()
    label = coh.map { c -> c.label }
    joint = coh.map { c -> c.joint }

    // Destructured: a record in a process `input:` is hashed by object identity, so its
    // task hash changes every run and -resume never matches it.
    lof = svannotate(label, coh.map { c -> c.vcf }, coh.map { c -> c.tbi },
                     protein_coding_gtf, noncoding_bed, symbolic_awk)
    sexed = sex_stratified_af(label, lof.out, ped)
    talos = talos_schema(label, sexed.out, schema_awk, keep_awk, check_awk, gnomad_pop, joint)

    all_versions = lof.versions.mix(sexed.versions, talos.versions)

    emit:
    vcf: Channel<MergedVcf> = talos.out
    report: Channel<Path> = talos.report
    versions: Channel<Path> = all_versions
}
