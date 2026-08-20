nextflow.enable.types = true
include { MergedVcf } from '../merge/svdb_merge.nf'

// Carries the sample name only, never its alignment. A Path rebuilt into an output record is
// a TaskPath naming a file in THAT task's work dir: feeding one back into a process throws
// UnsupportedOperationException out of the task hasher. The alignment therefore travels
// straight from the samplesheet channel to duphold, joined to these sites on the sample name
// -- see depth_svs.nf. Joining on the name is also stronger than the positional pairing this
// record used to provide.
record SampleSites {
    sample: String
    vcf: Path
    tbi: Path
}

// Phase 3: depth evidence.
//
// duphold ANNOTATES, it does not genotype. It attaches depth fold-change to calls that
// already exist -- DHFFC against 1 kb flanks, DHBFC against GC-matched bins -- which is
// what makes a ./. in a merged CNV row interpretable as real reference rather than a missed
// call. That ambiguity is load-bearing here because the design chose to keep multi-sample VCFs
// whole, and duphold is what earns the information back for CNVs.
//
// It is not a regenotyper. Paragraph and GraphTyper2 are, and both stay out of v1:
// stating the limitation is better than implying the genotypes are joint when they are not.

// Three processes rather than one, because each container only has its own tool -- the
// duphold image ships no bcftools, which is an exit 127 the moment you assume otherwise.
// Splitting is required regardless: duphold expects the VCF to hold exactly the one sample
// its alignment belongs to.

process split_sample {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 2.GB * task.attempt }

    input:
    cohort_vcf: Path
    cohort_tbi: Path
    sample: String
    // Lists SVDB's per-input blobs off the cohort header. Same asset the AnnotSV branch
    // uses, for the same reason in a different place -- see the script below.
    blob_awk: Path

    output:
    out: SampleSites = new SampleSites(
        sample: sample,
        vcf: file("${sample}.sites.vcf.gz"),
        tbi: file("${sample}.sites.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # SVDB's per-input blobs are stripped before the split, and this is a space decision, not
    # a cosmetic one. SVDB copies each input's whole source record into INFO as <tag>_INFO,
    # <tag>_SAMPLE and friends, on every row -- so a 300-input cohort carries ~300 source
    # records per record. Kept here, they ride through duphold into an UNCOMPRESSED
    # <sample>.dh.vcf, once per sample: tens of GB each, multi-TB of work dir for one cohort.
    #
    # Nothing can lose information by dropping them on this branch. duphold reads the
    # interval and SVTYPE; recombine_depth harvests only ID/SAMPLE/DH* back out of the
    # per-sample files and joins them onto the FULL cohort VCF, which keeps every key. The
    # per-input provenance stays answerable from the published VCF, exactly as it does for the
    # AnnotSV branch that strips the same keys for the same reason.
    #
    # Read off the header rather than named in a param: the key names carry each input's tag.
    # Empty on a cohort that never went through a merge, hence the guard.
    blobs=\$(bcftools view -h "${cohort_vcf.name}" | awk -f "${blob_awk.name}")
    if [ -n "\$blobs" ]; then
        bcftools view -s "${sample}" -Ou "${cohort_vcf.name}" \\
            | bcftools annotate -x "\$blobs" -Oz -o "${sample}.sites.vcf.gz"
    else
        bcftools view -s "${sample}" -Oz -o "${sample}.sites.vcf.gz" "${cohort_vcf.name}"
    fi
    bcftools index -t "${sample}.sites.vcf.gz"

    # Asserted, not assumed: a strip that quietly did nothing costs no correctness and gives
    # no error, it just puts the multi-TB work dir back. Cheap -- this reads the header only.
    if bcftools view -h "${sample}.sites.vcf.gz" | awk -f "${blob_awk.name}" | grep -q .; then
        echo "ERROR: SVDB blob INFO keys survived the strip in ${sample}.sites.vcf.gz" >&2
        exit 1
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// One sample, one alignment, one pass.
process duphold {
    container 'quay.io/biocontainers/duphold:0.2.1--h031d066_4'
    cpus 2
    memory { 8.GB * task.attempt }

    input:
    sample: String
    sites_vcf: Path
    sites_tbi: Path
    alignment: Path
    alignment_index: Path
    reference: Path
    reference_index: Path

    output:
    out: Path = file("${sample}.dh.vcf")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # These arrive as separate inputs rather than one record, because a record in an `input:`
    # block is hashed by object identity and defeats -resume. The sites VCF and the alignment
    # reach this process down two different routes and are joined on the sample name upstream
    # (depth_svs.nf); depth landing on the wrong sample is silent and severe, so the join is
    # checked rather than trusted. split_sample names the file after the sample it extracted.
    if [ "${sites_vcf.name}" != "${sample}.sites.vcf.gz" ]; then
        echo "ERROR: sites VCF ${sites_vcf.name} does not belong to sample ${sample}" >&2
        exit 1
    fi

    # CRAM cannot be decoded without the exact FASTA it was written against, and the wrong
    # one gives wrong depth rather than an error -- which then drives the depth soft filter. This is
    # params.alignment_reference, deliberately a separate knob from any annotation build.
    duphold \\
        --threads ${task.cpus} \\
        -v "${sites_vcf.name}" \\
        -b "${alignment.name}" \\
        -f "${reference.name}" \\
        -o "${sample}.dh.vcf"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    duphold: \$(duphold --help 2>&1 | head -n1 | sed 's/^ *version: //')
END_VERSIONS
    """
}

// Put the per-sample fields back, joined on record ID. Never positionally, never on
// POS/REF/ALT -- two symbolic-ALT records can share CHROM/POS/REF/ALT and differ only in
// INFO/END, and bcftools merge collapses those.
process recombine_depth {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    label: String
    vcf: Path
    tbi: Path
    joint: Boolean
    annotated: List<Path>
    extract_awk: Path
    merge_awk: Path

    output:
    out: MergedVcf = new MergedVcf(
        label: label,
        joint: joint,
        vcf: file("${label}.depth.vcf.gz"),
        tbi: file("${label}.depth.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    def dh_files = annotated.collect { t -> t.name }.join(' ')
    // A CONSTANT, deliberately not derived from task.memory. Nextflow leaves `memory` and
    // `cpus` out of the task hash precisely so resources can be tuned without invalidating
    // the cache -- but interpolating one INTO the script puts it back in, through the script
    // text. Sizing this from task.memory meant a `withName: recombine_depth { memory = ... }`
    // override re-ran the most expensive process in the pipeline, which is the opposite of
    // what reaching for that override is for.
    //
    // ponytail: fixed 1 GB buffer, not adaptive. These sorts are disk-bound -- a bigger
    // buffer buys roughly one fewer merge pass over the depth stream, not a different order
    // of magnitude -- and 1 GB sits safely inside this process's declared memory however it
    // is overridden. Raise this line if the sort phase ever dominates a real run; that will
    // rehash the task, which is correct, because it is a change to what the task does.
    def sort_mem = "1G"
    """
    set -euo pipefail
    # LC_ALL=C is load-bearing, not hygiene. merge_depth.awk is a merge join over two
    # ID-sorted streams, and gawk uses strcoll for `<` outside the C locale -- a collation
    # disagreeing with sort's byte order would misalign the streams and silently drop depth.
    # Both sorts and the awk run under the same collation here.
    export LC_ALL=C
    TAB=\$(printf '\\t')

    # The depth rows are piped straight into an ID-ordered stream rather than accumulated in
    # a file first. The unsorted intermediate this used to write reached 72 GB on a
    # 300-sample cohort, and the awk that read it into a hash needed an estimated 150 GB of
    # RAM -- it thrashed for 15 hours instead of failing. Nothing is materialised now beyond
    # the sorted stream and sort's own spill files.
    #
    # The extract runs here rather than in the duphold process because that image has no
    # bcftools; extract_depth.awk reads the FORMAT column by name so a tag duphold did not
    # emit becomes "." instead of a hard bcftools error.
    #
    # A cohort where no sample has an alignment is a supported configuration: the loop
    # then contributes nothing, the depth stream is empty, and merge_depth.awk pads every
    # column rather than failing.
    for f in ${dh_files}; do
        awk -f "${extract_awk.name}" "\$f"
    done | sort -t "\$TAB" -k1,1 -T . -S ${sort_mem} > depth.sorted.tsv

    # The join needs the VCF body in ID order too. Header and body are split so only the body
    # is sorted, then bcftools sort puts coordinate order back -- the same disk-based sort
    # promote_support and finalize already rely on.
    bcftools view -h "${vcf.name}" > header.vcf
    bcftools view -H "${vcf.name}" \\
        | sort -t "\$TAB" -k3,3 -T . -S ${sort_mem} > body.idsorted.vcf

    cat header.vcf body.idsorted.vcf \\
        | awk -f "${merge_awk.name}" depth.sorted.tsv - \\
        | bcftools sort -T . -Oz -o "${label}.depth.vcf.gz"
    bcftools index -t "${label}.depth.vcf.gz"

    # The join-correctness assertions, run rather than trusted: depth landing on the wrong record or the
    # wrong sample is silent and severe, so both are checked before anything is published.
    before=\$(bcftools view -H "${vcf.name}" | wc -l | tr -d ' ')
    after=\$(bcftools view -H "${label}.depth.vcf.gz" | wc -l | tr -d ' ')
    if [ "\$before" != "\$after" ]; then
        echo "ERROR: record count changed during depth recombination: \$before -> \$after" >&2
        exit 1
    fi
    bcftools query -l "${vcf.name}" > "before_samples.txt"
    bcftools query -l "${label}.depth.vcf.gz" > "after_samples.txt"
    if ! diff -q "before_samples.txt" "after_samples.txt" > /dev/null; then
        echo "ERROR: sample order changed during depth recombination" >&2
        exit 1
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}
