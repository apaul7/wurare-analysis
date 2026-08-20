nextflow.enable.types = true

// Phase 3b: per-sample ploidy, relatedness and (optionally) ancestry, from the alignments.
//
// WHY THIS EXISTS. The depth filter was autosome-only because a ploidy-aware threshold needs
// per-sample sex, and params.ped cannot be trusted to supply it -- a mislabelled PED row and
// a swapped sample look identical from inside the pipeline. somalier genotypes ~17k common
// sites through the alignment index and answers all three questions from one pass: X and Y
// depth (ploidy), pairwise relatedness (does the PED describe this cohort), and ancestry
// against a labelled reference.
//
// samtools idxstats was the obvious cheaper answer and was rejected: it is index-only on BAM
// but reads the ENTIRE file on CRAM ("does so by reading through the entire file" --
// samtools-idxstats(1)), it cannot filter on MAPQ, so chrX/chrY multi-mappers inflate a
// female's chrY signal, and it answers the ploidy question alone.
//
// Notably NOT `somalier relate --infer` when a pedigree is given: the ploidy call is derived
// from the depth columns in assets/check_somalier.awk instead. --infer changes which column
// carries the inferred sex, and depending on an inference flag's output shape is a worse
// dependency than three stable column names. --infer is also about inferring PEDIGREE, which
// is the operator's to state when they have stated it -- somalier is then only asked to check
// it, and inference must never compete with a stated pedigree. Only when the operator supplies
// NO pedigree does relate run with --infer: somalier's inference is then the sole source of
// family structure for the draft ped, while ploidy and sex are still derived from the depth
// columns in check_somalier.awk. One caveat: where depth alone cannot decide a sex,
// check_somalier keeps whatever the pedigree carried, so under --infer that residue is
// somalier's own het-ratio call rather than a depth measurement.

// One sample, one alignment. Cheap on both BAM and CRAM: somalier seeks to the sites rather
// than streaming the file, which is why this stage is affordable at cohort scale.
process somalier_extract {
    container 'quay.io/biocontainers/somalier:0.2.19--h0c29559_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    sample: String
    alignment: Path
    alignment_index: Path
    sites: Path
    reference: Path
    reference_index: Path

    output:
    out: Path = file("${sample}.somalier")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # CRAM cannot be decoded without the exact FASTA it was written against, and the wrong
    # one gives wrong depth rather than an error -- same knob, same hazard, as duphold.
    somalier extract \\
        -d . \\
        --sites "${sites.name}" \\
        -f "${reference.name}" \\
        "${alignment.name}"

    # somalier names its output after the @RG SM in the alignment, NOT after the file. A
    # sheet whose `sample` column disagrees with the read group is exactly the mislabelling
    # this stage exists to catch, and letting it through would attach one sample's ploidy to
    # another's calls -- silent, and severe.
    if [ ! -f "${sample}.somalier" ]; then
        echo "ERROR: somalier wrote \$(ls *.somalier 2>/dev/null | tr '\\n' ' ')rather than" \\
             "${sample}.somalier -- the alignment's read-group sample name does not match" \\
             "the alignments sheet for ${sample}" >&2
        exit 1
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    somalier: \$(somalier --version 2>&1 | head -n1 | sed 's/^somalier version: //')
END_VERSIONS
    """
}

// The same extraction from an already-merged multi-sample VCF, for a pipeline that is handed
// genotypes rather than alignments (annotate_snps). One task for the whole cohort, because
// somalier writes one .somalier per sample in the VCF.
//
// Depth here comes from FORMAT/AD rather than from reads, so this path is only as good as
// what the caller wrote: a VCF without AD yields zero depth, which check_somalier.awk turns
// into a hard error naming the cause rather than a cohort of undetermined samples.
// Relatedness needs only the genotypes and is unaffected.
process somalier_extract_vcf {
    container 'quay.io/biocontainers/somalier:0.2.19--h0c29559_0'
    cpus 1
    memory { 8.GB * task.attempt }

    input:
    vcf: Path
    tbi: Path
    sites: Path
    reference: Path
    reference_index: Path

    output:
    out: List<Path> = files("extracted/*.somalier")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    mkdir -p extracted
    somalier extract \\
        -d extracted \\
        --sites "${sites.name}" \\
        -f "${reference.name}" \\
        "${vcf.name}"

    # A VCF somalier could not genotype leaves an EMPTY directory rather than failing, and
    # `somalier relate` on nothing then dies with something unrelated to the cause.
    if ! ls extracted/*.somalier >/dev/null 2>&1; then
        echo "ERROR: somalier extracted nothing from ${vcf.name} -- the sites file and the" \\
             "VCF are most likely built against different references" >&2
        exit 1
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    somalier: \$(somalier --version 2>&1 | head -n1 | sed 's/^somalier version: //')
END_VERSIONS
    """
}

// One task over the whole cohort. Relatedness is pairwise, so it cannot be sharded.
process somalier_relate {
    container 'quay.io/biocontainers/somalier:0.2.19--h0c29559_0'
    cpus 1
    memory { 8.GB * task.attempt }

    input:
    label: String
    extracted: List<Path>
    // Empty list, or exactly one PED. Optional because annotate_snps' QC is reporting-only:
    // a cohort with no pedigree still gets useful observed relatedness, and demanding a file
    // whose only job is to be contradicted would block the check entirely. annotate_svs always
    // passes one -- it consumes the inferred PED downstream.
    // Not a Path: typed syntax has no null, and an empty list is how "absent" is expressed for
    // an optional file input -- the same idiom as annotate_snps' precomputed_cadd.
    ped: List<Path>
    // Only meaningful when ped is empty: asks somalier to infer family structure so the
    // draft ped below carries something better than singletons. Never set when a PED is
    // given -- inference must not compete with a stated pedigree.
    infer: Boolean

    output:
    samples: Path = file("${label}.samples.tsv")
    pairs: Path = file("${label}.pairs.tsv")
    html: Path = file("${label}.html")
    draft_ped: Path = file("${label}.draft.ped")
    versions: Path = file("versions.yml")

    script:
    def files = extracted.collect { t -> t.name }.join(' ')
    def ped_arg = ped ? "--ped ${ped[0].name}" : ""
    def infer_arg = infer ? "--infer" : ""
    """
    set -euo pipefail
    # --ped, when given, so pairs.tsv carries expected_relatedness and the observed values can
    # be checked against what the pedigree claims. Without one, --infer asks somalier to infer
    # that structure instead, and the draft ped below carries it onward. Sex is NOT taken from
    # here either way -- see the header comment.
    somalier relate ${infer_arg} ${ped_arg} -o "${label}" ${files}

    # samples.tsv's first six columns are ped-shaped (#family_id..phenotype). Cut them out,
    # header stripped, as the draft pedigree check_somalier consumes when no operator PED
    # exists. Built unconditionally: typed outputs have no optional, and an unused file is
    # cheaper than one. somalier writes -9 for unknown parent/sex where the PED convention
    # downstream tools expect is 0, so those columns are normalized; phenotype keeps -9.
    awk -F'\\t' -v OFS='\\t' '!/^#/ {
        for (i = 3; i <= 5; i++) if (\$i == "-9") \$i = 0
        print \$1, \$2, \$3, \$4, \$5, \$6
    }' "${label}.samples.tsv" > "${label}.draft.ped"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    somalier: \$(somalier --version 2>&1 | head -n1 | sed 's/^somalier version: //')
END_VERSIONS
    """
}

// Optional, and separate because it needs a labelled reference cohort the other two do not.
// QC in its own right: a predicted ancestry contradicting the referral is the same class of
// signal as a relatedness mismatch, and it is what makes a population-AF threshold defensible
// for a non-European case.
process somalier_ancestry {
    container 'quay.io/biocontainers/somalier:0.2.19--h0c29559_0'
    cpus 1
    memory { 8.GB * task.attempt }

    input:
    label: String
    extracted: List<Path>
    labels: Path
    labelled_dir: Path

    output:
    // Every row somalier wrote, the labelled reference samples included. That is what makes
    // the PCA readable -- a predicted ancestry means little without the background it was
    // predicted against.
    out: Path = file("somalier-ancestry.tsv")
    // The same table narrowed to this run's samples. The reference rows outnumber a cohort by
    // orders of magnitude, so the full table is the wrong thing to read a result out of.
    cohort_out: Path = file("${label}.somalier-ancestry.tsv")
    html: Path = file("somalier-ancestry.html")
    versions: Path = file("versions.yml")

    script:
    def files = extracted.collect { t -> t.name }.join(' ')
    // The run's sample ids, taken from the staged filenames -- somalier extract writes
    // <sample>.somalier, and its #sample_id column carries that same sample name.
    def cohort_ids = extracted.collect { t -> t.name.replaceFirst(/\.somalier$/, '') }.join(',')
    """
    set -euo pipefail
    # The labelled reference set goes before ++, this cohort after it. Reversing them
    # predicts the reference cohort's ancestry from ours, which is not an error anywhere.
    somalier ancestry \\
        --labels "${labels.name}" \\
        "${labelled_dir.name}"/*.somalier \\
        ++ ${files}

    # somalier appends a FIXED ".somalier-ancestry.<ext>" to its output prefix, and that prefix
    # defaults to "somalier-ancestry" -- so the files it actually writes are
    # somalier-ancestry.somalier-ancestry.tsv/.html. Renamed rather than given a different -o:
    # -o moves the prefix but never removes the suffix, so the undoubled name is only reachable
    # this way. Pinned in tests/test_somalier_assumptions.py against the real tool.
    mv somalier-ancestry.somalier-ancestry.tsv somalier-ancestry.tsv
    mv somalier-ancestry.somalier-ancestry.html somalier-ancestry.html

    # Filtered on this run's sample ids rather than on given_ancestry being empty. Both would
    # work today -- somalier leaves that column blank for query samples -- but the ids are
    # already in hand, where the column is somalier's to rename.
    awk -F'\\t' -v OFS='\\t' -v ids="${cohort_ids}" '
        BEGIN { n = split(ids, a, ","); for (i = 1; i <= n; i++) keep[a[i]] = 1 }
        NR == 1 { print; next }
        \$1 in keep { print }
    ' somalier-ancestry.tsv > "${label}.somalier-ancestry.tsv"

    # A cohort table with no rows means the .somalier filenames and somalier's #sample_id
    # column disagree, which would otherwise publish an empty file and look like a result.
    if [ "\$(awk 'NR > 1' "${label}.somalier-ancestry.tsv" | wc -l)" -eq 0 ]; then
        echo "ERROR: none of this run's samples (${cohort_ids}) appear in somalier-ancestry.tsv" >&2
        exit 1
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    somalier: \$(somalier --version 2>&1 | head -n1 | sed 's/^somalier version: //')
END_VERSIONS
    """
}

// somalier's tables -> the ploidy table the depth filter reads, plus a PED carrying the
// inferred sex. The awk holds the thresholds and every assertion; see
// assets/check_somalier.awk.
process check_somalier {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 2.GB * task.attempt }

    input:
    samples: Path
    pairs: Path
    ped: Path
    check_awk: Path

    output:
    ploidy: Path = file("ploidy.tsv")
    inferred_ped: Path = file("cohort.inferred.ped")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    awk -v SAMPLES="${samples.name}" \\
        -v PAIRS="${pairs.name}" \\
        -v PED="${ped.name}" \\
        -v OUT_PLOIDY=ploidy.tsv \\
        -v OUT_PED=cohort.inferred.ped \\
        -f "${check_awk.name}" </dev/null

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    awk: \$(awk --version 2>&1 | head -n1)
END_VERSIONS
    """
}
