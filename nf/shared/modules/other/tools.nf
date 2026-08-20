// Copied from the upstream pipeline this was ported from.
//
// NO LONGER VERBATIM, beyond the repo-wide container tag pin: split_vcf's `input:` takes
// scalars, not a record, and the `stub:` blocks are absent (this repo tests by running the
// real DAG against fake tools instead -- see tests/run.sh).
// Everything else is unchanged. The scalar `input:` is upstream's defect too.
nextflow.enable.types = true

// Lists the distinct chromosomes actually present in a VCF, so annotate_variants can
// split/annotate only regions with real data instead of a static chr1..chrY list --
// on targeted-panel data most of those chromosomes have zero variants, and CADD/
// SpliceAI/SQUIRLS/VEP all error out when handed an empty VCF.
process list_chroms {
    container "apaul7/analysis:1.2.0"
    cpus 1
    memory { 1.GB * task.attempt }

    input:
    tuple(vcf: Path, tbi: Path)

    output:
    out:      Path = file("chroms.txt")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    bcftools query --format '%CHROM\\n' "${vcf.name}" | sort -u > chroms.txt

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}
// Filename-safe token for an interval string. Contig names are not guaranteed to be safe
// to put in a filename: ALT/HLA contigs carry `*` and `:` (e.g. HLA-A*01:01:01:01), and a
// region string like chr1:1000-2000 has the same problem. Anything outside
// [A-Za-z0-9._-] becomes an underscore.
//
// This is deliberately NOT reversible, which is exactly why `interval` and `tag` are
// carried as two separate fields from split_vcf onward: `interval` stays the raw string
// that bcftools -r and pysam's fetch() need, while `tag` is used for every output
// filename. Deriving one from the other -- or recovering either from `vcf.simpleName`,
// which also truncates at the first dot -- is what this split exists to prevent.
def interval_tag(interval) {
    return interval.replaceAll(/[^A-Za-z0-9._-]/, '_')
}

process split_vcf {
    container "mgibio/bcftools-cwl:1.12"
    cpus 1
    memory { 10.GB * task.attempt }

    // Scalars, not a record. A record in `input:` has no value-based hash, so the task hasher
    // falls back to Object.toString() and embeds the JVM identity hash -- a fresh address
    // every run, so `-resume` never matches. Only `input:` is affected; the `output:` record
    // below is fine. Canonical statement for annotate_snps; annotate_svs has its own in
    // modules/merge/svdb_merge.nf.
    input:
    tuple(vcf: Path, tbi: Path)
    interval: String
    tag: String

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        vcf: file("${tag}.out.vcf.gz"),
        tbi: file("${tag}.out.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # -Oz is required: bcftools does not infer the output format from the .gz extension.
    bcftools view -r "${interval}" -Oz -o "${tag}.out.vcf.gz" "${vcf.name}"
    bcftools index -t "${tag}.out.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}
