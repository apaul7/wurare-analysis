nextflow.enable.types = true

// One tidy TSV (#section  metric  sample  value) summarizing the whole run -- variant
// counts, consequence breakdowns, inferred sex/ancestry -- for a cursory sanity glance
// without opening the multi-GB multianno TSV or the somalier tables.
//
// Stats run against the INPUT VCF, not the annotated one: the records are identical
// (SpliceAI/SQUIRLS only add INFO), and the input exists even under
// --skip_spliceai_squirls, so the report needs no fork on that flag.
//
// The multianno TSV is the RARE subset only -- reduce_variants filters to gnomAD
// AAF < 0.01 before table_annovar -- which is why the report's whole-cohort numbers must
// come from bcftools stats and the consequence counts are labeled rare_subset.
//
// summary_report.py is deliberately tolerant (missing somalier files, renamed columns ->
// rows omitted with a stderr note, exit 0): a report must never kill the run that
// produced the data it summarizes.
process summary_report {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 5.GB * task.attempt }

    input:
    vcf: Path
    tbi: Path
    multianno: Path
    // Empty list, or exactly one file each -- the QC stage is opt-in and the PED optional,
    // so any of these can be absent. Same idiom as somalier_relate's ped: typed syntax has
    // no null, and an empty list is how "absent" is expressed for an optional file input.
    samples_tsv: List<Path>
    ploidy: List<Path>
    ancestry: List<Path>
    pairs: List<Path>
    cohort: String
    data_type: String
    run_date: String
    clinvar_date: String

    output:
    out: Path = file("${cohort}_${data_type}_${run_date}.summary.tsv")
    versions: Path = file("versions.yml")

    script:
    def samples_arg  = samples_tsv ? "--samples-tsv \"${samples_tsv[0].name}\"" : ""
    def ploidy_arg   = ploidy ? "--ploidy \"${ploidy[0].name}\"" : ""
    def ancestry_arg = ancestry ? "--ancestry \"${ancestry[0].name}\"" : ""
    def pairs_arg    = pairs ? "--pairs \"${pairs[0].name}\"" : ""
    """
    set -euo pipefail
    bcftools stats -s - "${vcf.name}" > bcftools_stats.txt

    summary_report.py \\
        --stats bcftools_stats.txt \\
        --multianno "${multianno.name}" \\
        --cohort "${cohort}" \\
        --data-type "${data_type}" \\
        --run-date "${run_date}" \\
        --clinvar-date "${clinvar_date}" \\
        ${samples_arg} \\
        ${ploidy_arg} \\
        ${ancestry_arg} \\
        ${pairs_arg} \\
        --out "${cohort}_${data_type}_${run_date}.summary.tsv"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    python: \$(python3 --version 2>&1 | sed 's/^Python //')
END_VERSIONS
    """
}
