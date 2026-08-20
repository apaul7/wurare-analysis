nextflow.enable.types = true

// Path inputs are staged into the task work dir, so scripts refer to them by `.name` -- the
// absolute host path is not visible inside the container. Same convention as annotate_snps'
// squirls.nf.
//
// These records carry the sheet row's three identifying scalars flat, rather than nesting
// the whole SvInput. They used to nest it, and the nested `vcf`/`tbi` were unusable: a Path
// rebuilt into an output record is a TaskPath naming a file in THAT task's work dir, so
// feeding one back into a process throws UnsupportedOperationException out of the task
// hasher. Nothing ever read them -- every consumer wants sample_set, caller and joint -- so
// they are gone rather than carried broken.
record SampleIds {
    sample_set: String
    caller: String
    joint: Boolean
    samples: Path
}

record FilteredVcf {
    sample_set: String
    caller: String
    joint: Boolean
    vcf: Path
}

record StandardVcf {
    sample_set: String
    caller: String
    joint: Boolean
    vcf: Path
    tbi: Path
}

// The standardized VCF plus the key the merge groups on. That key is the sorted set of
// sample IDs read from the VCF header -- never the sheet's sample_set label, which is free
// text a human types and by design is NOT the grouping key. Carried as its own
// record rather than recomputed downstream so there is exactly one definition of it.
record GroupedVcf {
    sample_key: String
    entry: StandardVcf
}

// Phase 1: get every input VCF onto one schema before any merge sees it --
// merging raw caller output is where silent corruption enters.
//
// Kept in one file rather than the three the design sketched (size_filter/standardize/stamp_ids):
// they are one linear stage over two containers, and splitting them buys three include
// lines and no clarity. Split when one of them grows a reason to vary on its own.

// Sample IDs come from the VCF header, never from the sheet, so the derived grouping
// and the join-key check are based on what the files actually contain.
process read_sample_ids {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 2.GB * task.attempt }

    input:
    sample_set: String
    caller: String
    joint: Boolean
    vcf: Path
    tbi: Path

    output:
    out: SampleIds = new SampleIds(
        sample_set: sample_set,
        caller: caller,
        joint: joint,
        samples: file("${sample_set}.${caller}.samples.txt")
    )
    versions: Path = file("versions.yml")

    script:
    def tag = "${sample_set}.${caller}"
    """
    set -euo pipefail
    bcftools query -l "${vcf.name}" > "${tag}.samples.txt"
    if [ ! -s "${tag}.samples.txt" ]; then
        echo "ERROR: ${vcf.name} has no sample columns -- a sites-only VCF cannot" >&2
        echo "       contribute genotypes to the cohort matrix" >&2
        exit 1
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// PASS,. and not PASS: a bare -f PASS drops every record whose FILTER is '.', which is how
// several callers spell "unfiltered". That is a total, silent loss of those inputs, and it
// surfaces three stages later looking like a merge problem. Counts are recorded either
// side of the filter so a caller that loses everything shows up as a round zero rather than
// as an absence nobody notices.
//
// BND is exempt from the size floor: a breakend has no meaningful SVLEN, so filtering it on
// one drops translocations silently -- the classic BND failure mode.
process pass_size_filter {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    sample_set: String
    caller: String
    joint: Boolean
    vcf: Path
    tbi: Path
    min_sv_size: Integer

    output:
    out: FilteredVcf = new FilteredVcf(
        sample_set: sample_set,
        caller: caller,
        joint: joint,
        vcf: file("${sample_set}.${caller}.filtered.vcf.gz")
    )
    counts: Path = file("${sample_set}.${caller}.filter_counts.tsv")
    versions: Path = file("versions.yml")

    script:
    def tag = "${sample_set}.${caller}"
    """
    set -euo pipefail
    before=\$(bcftools view -H "${vcf.name}" | wc -l | tr -d ' ')

    bcftools view -f PASS,. -Ou "${vcf.name}" \\
        | bcftools filter \\
            -e 'INFO/SVTYPE!="BND" && INFO/SVLEN!="." && abs(INFO/SVLEN) < ${min_sv_size}' \\
            -Oz -o "${tag}.filtered.vcf.gz"

    after=\$(bcftools view -H "${tag}.filtered.vcf.gz" | wc -l | tr -d ' ')
    printf 'sample_set\\tcaller\\tbefore\\tafter\\n' > "${tag}.filter_counts.tsv"
    printf '%s\\t%s\\t%s\\t%s\\n' "${sample_set}" "${caller}" "\$before" "\$after" \\
        >> "${tag}.filter_counts.tsv"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// Normalize, sort, stamp a unique record ID, stamp ALGORITHMS, compress and index.
//
// Normalization is assets/normalize_records.awk and applies to every input, with no
// caller branch. `svtk standardize` used to own this and is gone: it cannot run at all
// (invalid bundled template under modern pysam, every biocontainer build), and keeping a
// branch for a tool that never executes is a branch that only ever misleads.
//
// The ID is load-bearing later, not cosmetic: duphold's per-sample FORMAT recombination
// rejoins per-sample outputs onto the same records, and positional or POS/REF/ALT
// joins on symbolic SV alleles are exactly how depth lands on the wrong record. It is also
// what `bcftools merge -m id` needs for the merge's truvari upgrade path. A running counter rather
// than CHROM_POS, because two callers can and do emit two records at one position.
//
// The stamping itself is assets/stamp_records.awk, staged as an input rather than written
// inline. An awk program inside a Nextflow triple-quoted string has to survive Groovy
// interpolation and then the shell, and a stray dollar or quote in that string silently
// shifts every later interpolation by one position -- which is exactly how the first two
// attempts at this broke. As a staged file it is plain awk, and
// tests/test_prepare_stamps.py exercises it with no Nextflow and no container.
//
// ALGORITHMS is stamped here because Talos defaults it to gCNV when absent, so an unstamped
// callset is mislabelled rather than rejected. It is re-derived after each merge --
// SVDB keeps only the priority record's value -- so this stamp is the input to that
// derivation, never the final answer.
process finalize {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    sample_set: String
    caller: String
    joint: Boolean
    filtered_vcf: Path
    normalize_awk: Path
    stamp_awk: Path

    output:
    out: StandardVcf = new StandardVcf(
        sample_set: sample_set,
        caller: caller,
        joint: joint,
        vcf: file("${sample_set}.${caller}.std.vcf.gz"),
        tbi: file("${sample_set}.${caller}.std.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    def tag = "${sample_set}.${caller}"
    """
    set -euo pipefail
    # Staged from assets/ rather than written inline -- see the note above the process.
    bcftools sort -Ov "${filtered_vcf.name}" \\
        | awk -f "${normalize_awk.name}" \\
        | awk -v TAG="${tag}" -v CALLER="${caller}" -f "${stamp_awk.name}" \\
        | bgzip -c > "${tag}.std.vcf.gz"

    bcftools index -t "${tag}.std.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}
