nextflow.enable.types = true
include { MergedVcf } from '../merge/svdb_merge.nf'

// Phase 5: soft filter tags.
//
// Annotate hard, filter softly. Every record survives; failing a criterion only writes a
// tag into FILTER. One run then serves Talos, manual review and a QC report without
// re-running the expensive half, and each consumer picks its own threshold from the tags
// rather than from a filtered file it cannot widen again.
//
// The thresholds live in assets/tag_filters.awk, staged rather than inline -- the same
// reason as everywhere else in this pipeline: a stray dollar or quote in a Nextflow
// triple-quoted script shifts every later interpolation by one position, silently.
process tag_filters {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    label: String
    vcf: Path
    tbi: Path
    joint: Boolean
    filter_awk: Path
    // Per-sample chrX/chrY copy numbers, from somalier (nf/shared/assets/check_somalier.awk). Always a
    // real file: with no somalier stage it is assets/ploidy_none.tsv, which names no sample,
    // and the awk then leaves the sex chromosomes exempt exactly as it did before.
    ploidy: Path
    // Only ever used to pick PAR coordinates. An unrecognized build falls back to exempting
    // chrX/chrY rather than guessing at them.
    genome_build: String
    pop_af: String
    internal_af: String
    del_dhffc: String
    dup_dhbfc: String
    min_callers: String

    output:
    out: MergedVcf = new MergedVcf(
        label: label,
        joint: joint,
        vcf: file("${label}.tagged.vcf.gz"),
        tbi: file("${label}.tagged.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    bcftools view "${vcf.name}" \\
        | awk -v POP_AF="${pop_af}" -v INT_AF="${internal_af}" \\
              -v DEL_DHFFC="${del_dhffc}" -v DUP_DHBFC="${dup_dhbfc}" \\
              -v MIN_CALLERS="${min_callers}" \\
              -v PLOIDY="${ploidy.name}" -v BUILD="${genome_build}" \\
              -f "${filter_awk.name}" \\
        | bcftools view -Oz -o "${label}.tagged.vcf.gz"
    bcftools index -t "${label}.tagged.vcf.gz"

    # Soft means soft: assert no record was lost. The contract is that nothing is removed from
    # the annotated VCF, only tagged, and an assertion is cheaper than trusting it.
    before=\$(bcftools view -H "${vcf.name}" | wc -l | tr -d ' ')
    after=\$(bcftools view -H "${label}.tagged.vcf.gz" | wc -l | tr -d ' ')
    if [ "\$before" != "\$after" ]; then
        echo "ERROR: filtering removed records (\$before -> \$after); it must only tag" >&2
        exit 1
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}
