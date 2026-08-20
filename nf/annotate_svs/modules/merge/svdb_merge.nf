nextflow.enable.types = true

// Phase 2: the two-axis merge.
//
// One svdb_merge process, invoked at both axes -- axis A reconciles callers over the same
// sample set, axis B assembles the cohort. They differ only in their inputs, their tags and
// their thresholds, so a second process would be a copy that drifts.
//
// Path inputs are staged into the task work dir, so scripts use `.name`, never the absolute
// host path.
//
// Records are built and read in the workflows; they never appear in a process `input:`
// block. Nextflow has no value-based hash for a record, so it falls back to Object.toString
// -- which embeds the JVM identity hash. The task hash then changes on every run and
// -resume can never match the cache entry. Processes take label/vcf/tbi/joint instead, and
// the workflow destructures at the call site. Records in `output:` are unaffected and stay.
//
// And note: a bare dollar sign anywhere inside a
// triple-quoted script -- including in a comment -- starts a Groovy interpolation and
// silently shifts every later interpolation by one position. That cost two debugging rounds
// in Phase 1; it is why the INFO rewrites live in assets/*.awk instead of inline here.

// `joint` rides on the record because axis B's --priority order is decided by it: a
// sample present in both a joint VCF and its own single-sample calls keeps the joint
// genotype only if the joint input is named first, and that was measured, not assumed
// (in the design spike). Carried as a field rather than rejoined by label downstream -- a second channel
// keyed on a filename is the plumbing that goes quietly wrong. Past axis B it is inert:
// every stage from depth onward copies it through, because the cohort VCF is joint by then
// and nothing reads the flag again.
record MergedVcf {
    label: String
    joint: Boolean
    vcf: Path
    tbi: Path
}

// svdb_merge's output carries its own label rather than emitting a bare Path.
//
// Not cosmetic. A process output channel emits in TASK COMPLETION order, so pairing it with a
// label taken from a separate channel -- which emits in channel order -- glues the nth label
// to the nth-FINISHED file. promote_axis_a did exactly that, and the consequences were two:
// its task hash changed on every run so -resume never matched it, and worse, the pairing is a
// permutation, so a group's merged VCF could carry another group's label and another group's
// `joint`. The label becomes the axis-B tag and `joint` decides axis-B priority, which the
// design spike measured as changing which genotype survives. With one group it cannot manifest, which is
// why the fixtures never caught it.
//
// Carrying the label on the item makes the downstream join keyed rather than positional. Only
// `output:` blocks may hold records -- a record in an `input:` block is hashed by object
// identity and defeats -resume, which is why the processes below take primitives.
record AxisMerge {
    label: String
    vcf: Path
}

// SVDB needs `--vcf file:tag` notation for `--priority` to apply at all; a tag/priority
// mismatch is a hard error rather than a silent no-op (verified in the design spike), so the lists are
// built from one source and SVDB is left to complain if they ever drift.
//
// Inputs arrive already sorted by the caller of this process. Nextflow emits collected
// items in completion order, so without an explicit sort two runs over identical inputs can
// hand SVDB its files in different orders and produce different cohort VCFs.
process svdb_merge {
    container 'quay.io/biocontainers/svdb:2.12.0--py312hfcd9dac_0'
    cpus 1
    // Doubling, not `8.GB * task.attempt`. The two axes are one process with wildly different
    // appetites: axis A merges a single sample's callers and fits in 8 GB on the first try,
    // while axis B holds the whole cohort at once and grows with variants x samples. Linear
    // growth spent all three retries between 8 and 24 GB, so a cohort needing 40 GB failed
    // four times without ever asking for it. Doubling reaches 64 GB over the same retries and
    // costs a small cohort nothing, because it never leaves the first attempt.
    //
    // ponytail: headroom, not a fix for the shape. svdb merge is a single non-parallel process
    // over the entire cohort, so this climbs with N until it doesn't fit. Per-contig
    // sharding is the real answer and needs explicit BND handling first. Override per run
    // without editing this: process { withName: 'merge_svs:merge_axis_b' { memory = 128.GB } }
    memory { 8.GB * (2 ** (task.attempt - 1)) }

    input:
    label: String
    vcfs: List<Path>
    tags: List<String>
    overlap: String
    bnd_distance: String
    ins_distance: String

    output:
    out: AxisMerge = new AxisMerge(label: label, vcf: file("${label}.merged.vcf"))
    versions: Path = file("versions.yml")

    script:
    // Tags and priority come from one list, in one order, so they cannot disagree.
    def tagged = [vcfs, tags].transpose().collect { v, t -> "${v.name}:${t}" }.join(' ')
    def priority = tags.join(',')
    """
    set -euo pipefail
    svdb --merge \\
        --vcf ${tagged} \\
        --priority ${priority} \\
        --overlap ${overlap} \\
        --bnd_distance ${bnd_distance} \\
        --ins_distance ${ins_distance} \\
        > "${label}.merged.vcf"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    svdb: \$(svdb 2>&1 | head -n1 | sed 's/^usage: SVDB-//; s/:.*//')
END_VERSIONS
    """
}

// Mandatory between the axes, not cosmetic -- verified in the design spike. Also re-derives ALGORITHMS from the
// support, since SVDB keeps only the priority record's value.
process promote_support {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    label: String
    vcf: Path
    promote_awk: Path
    mode: String
    joint: Boolean

    output:
    out: MergedVcf = new MergedVcf(
        label: label,
        joint: joint,
        vcf: file("${label}.vcf.gz"),
        tbi: file("${label}.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # The label and the VCF arrive as separate inputs, so a rewiring that pairs them by
    # arrival order rather than by key mislabels a whole group -- silently, because every
    # group's data is still present, just permuted. The label becomes the axis-B tag and
    # travels with `joint`, which decides axis-B priority and therefore which genotype
    # survives. svdb_merge names its output after the label it merged, so the pairing
    # is checkable rather than assumed. This is what the zip bug looked like from here.
    if [ "${vcf.name}" != "${label}.merged.vcf" ]; then
        echo "ERROR: ${vcf.name} does not belong to label ${label}" >&2
        exit 1
    fi

    awk -v MODE="${mode}" -f "${promote_awk.name}" "${vcf.name}" "${vcf.name}" \\
        | bcftools sort -Oz -o "${label}.vcf.gz"
    bcftools index -t "${label}.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// Internal cohort AF, in the main line rather than the Talos tail. The soft filters use it
// and the design carries two risks about it, but before this stage nothing computed it -- compute_AFs.py
// lives in Phase 6, which most consumers never reach.
//
// The resulting AF has a hard floor at small N (a het singleton in 10 samples is AF 0.05)
// and is taken over a cohort of mixed genotype provenance. Both belong in the README, not
// in a comment nobody reads at review time.
process fill_tags {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 8.GB * task.attempt }

    input:
    label: String
    vcf: Path
    tbi: Path
    joint: Boolean

    output:
    out: MergedVcf = new MergedVcf(
        label: label,
        joint: joint,
        vcf: file("${label}.af.vcf.gz"),
        tbi: file("${label}.af.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    bcftools +fill-tags "${vcf.name}" -Oz -o "${label}.af.vcf.gz" -- -t AC,AN,AF
    bcftools index -t "${label}.af.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// A sample set with only one input has nothing to reconcile at axis A, but still has to
// reach axis B as an indexed bgzipped VCF like everything else. Repacking rather than
// passing the Phase 1 file straight through keeps axis B's inputs uniform.
process passthrough {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 2.GB * task.attempt }

    input:
    label: String
    vcf: Path
    joint: Boolean

    output:
    out: MergedVcf = new MergedVcf(
        label: label,
        joint: joint,
        vcf: file("${label}.vcf.gz"),
        tbi: file("${label}.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    bcftools view -Oz -o "${label}.vcf.gz" "${vcf.name}"
    bcftools index -t "${label}.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}
