// Copied from an upstream pipeline and locally modified; upstream is not tracked.
//
// NO LONGER VERBATIM, beyond the repo-wide container tag pin:
//   1. All three processes' `input:` blocks take scalars, not records -- see
//      nf/shared/modules/other/tools.nf. Upstream's defect too, as is (2).
//   2. The version is read from the module rather than hardcoded to a literal.
//   3. No `stub:` blocks (this repo tests by running the real DAG against fake tools).
//
// CONVERGED with upstream, not diverged: add_precomputed calls ../../bin/add_scores.py rather
// than heredoc'ing it, which is upstream's change. That script itself does carry
// one deliberate divergence, in how it treats a multi-allelic record -- read its header before
// syncing it either way.
nextflow.enable.types = true
include { Reference } from "../../../shared/modules/types.nf"
process run_spliceai {
    //label "gpu"
    container "apaul7/docker-splice-ai:v1.3.1.gpu"
    cpus 4
    memory { 50.GB * task.attempt }

    input:
    interval: String
    tag: String
    scored_vcf: Path
    scored_tbi: Path
    unscored_vcf: Path
    unscored_tbi: Path
    ref: Reference

    // The already-scored half is re-declared as an output so it stays paired with the
    // newly-scored half through to concat, rather than being re-joined by a second channel.
    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        newly_scored_vcf: file("${tag}.vcf"),
        prescored_vcf: file("${scored_vcf.name}"),
        prescored_tbi: file("${scored_tbi.name}")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    spliceai \\
        -A grch38 \\
        -R "${ref.fa}" \\
        -I "${unscored_vcf.name}" \\
        -O "${tag}.vcf"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    spliceai: \$(python3 -c 'import spliceai; print(spliceai.__version__)')
    python: \$(python3 --version 2>&1 | sed 's/^Python //')
END_VERSIONS
    """
}

process add_precomputed {
    container "apaul7/analysis:1.2.0"
    cpus 1
    memory { 10.GB * task.attempt }

    input:
    interval: String
    tag: String
    vcf: Path
    tbi: Path
    precomputed: Path
    // The index is a real input so it is staged next to `precomputed` and pysam can seek.
    // Previously this process ran `tabix` on the precomputed file itself, re-indexing a
    // genome-wide VCF inside every per-interval task.
    precomputed_tbi: Path

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        scored_vcf: file("${tag}.scored.vcf.gz"),
        scored_tbi: file("${tag}.scored.vcf.gz.tbi"),
        unscored_vcf: file("${tag}.unscored.vcf.gz"),
        unscored_tbi: file("${tag}.unscored.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    // add_scores.py lives in ../../bin/, which Nextflow puts on PATH automatically. It used
    // to be a ~60-line heredoc written to disk at runtime, which meant no syntax checking, no
    // linting and nothing that could test it in isolation.
    script:
    """
    set -euo pipefail
    add_scores.py \\
        --input "${vcf.name}" \\
        --precomputed "${precomputed}" \\
        --region "${interval}" \\
        --scored "${tag}.scored.vcf.gz" \\
        --unscored "${tag}.unscored.vcf.gz"
    tabix -p vcf "${tag}.scored.vcf.gz"
    tabix -p vcf "${tag}.unscored.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/^Python //')
    pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^tabix (htslib) //')
END_VERSIONS
    """
}
process concat {
    container "apaul7/analysis:1.2.0"
    cpus 1
    memory { 10.GB * task.attempt }

    input:
    interval: String
    tag: String
    newly_scored_vcf: Path
    prescored_vcf: Path
    prescored_tbi: Path

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        vcf: file("${tag}.spliceai.vcf.gz"),
        tbi: file("${tag}.spliceai.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    bgzip "${newly_scored_vcf.name}"
    tabix -p vcf "${newly_scored_vcf.name}.gz"

    # -Oz is required: bcftools does not infer the output format from the .gz extension.
    bcftools concat --allow-overlaps "${newly_scored_vcf.name}.gz" "${prescored_vcf.name}" | bcftools sort -Oz -o "${tag}.spliceai.vcf.gz"
    tabix -p vcf "${tag}.spliceai.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^tabix (htslib) //')
END_VERSIONS
    """
}
workflow spliceai {
    take:
    vcfs: Channel<Record>
    ref: Reference
    spliceai_precomputed_scores: Path
    spliceai_precomputed_tbi: Path

    main:
    // multiMap: one traversal, fields stay visibly one row -- see prepare_svs.nf.
    in_parts = vcfs.multiMap { r ->
        interval: r.interval
        tag: r.tag
        vcf: r.vcf
        tbi: r.tbi
    }
    add_precomputed_res = add_precomputed(in_parts.interval, in_parts.tag,
                                          in_parts.vcf, in_parts.tbi,
                                          spliceai_precomputed_scores, spliceai_precomputed_tbi)

    pre_parts = add_precomputed_res.out.multiMap { r ->
        interval: r.interval
        tag: r.tag
        scored_vcf: r.scored_vcf
        scored_tbi: r.scored_tbi
        unscored_vcf: r.unscored_vcf
        unscored_tbi: r.unscored_tbi
    }
    run_spliceai_res = run_spliceai(pre_parts.interval, pre_parts.tag,
                                    pre_parts.scored_vcf, pre_parts.scored_tbi,
                                    pre_parts.unscored_vcf, pre_parts.unscored_tbi, ref)

    scored_parts = run_spliceai_res.out.multiMap { r ->
        interval: r.interval
        tag: r.tag
        newly_scored_vcf: r.newly_scored_vcf
        prescored_vcf: r.prescored_vcf
        prescored_tbi: r.prescored_tbi
    }
    concat_res = concat(scored_parts.interval, scored_parts.tag,
                        scored_parts.newly_scored_vcf,
                        scored_parts.prescored_vcf, scored_parts.prescored_tbi)

    emit:
    vcf: Channel<Record> = concat_res.out
    versions: Channel<Path> = add_precomputed_res.versions
        .mix(run_spliceai_res.versions, concat_res.versions)
}
