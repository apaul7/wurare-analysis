// Copied from an upstream pipeline and locally modified; upstream is not tracked.
//
// NO LONGER VERBATIM, beyond the repo-wide container tag pin. Local hunks, all upstream
// defects too:
//   1. Both processes' `input:` blocks take scalars, not records -- see
//      nf/shared/modules/other/tools.nf.
//   2. The version is read from the jar rather than hardcoded to a literal.
// Also: no `stub:` blocks (this repo tests by running the real DAG against fake tools).
nextflow.enable.types = true
process run_squirls {
    container 'apaul7/docker-squirls:1.0.0'
    cpus 1
    memory { 25.GB * task.attempt }

    input:
    interval: String
    tag: String
    vcf: Path
    tbi: Path
    config: Path
    jannovar: Path

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        vcf: file("${tag}.spliceai.squirls.vcf")
    )
    versions: Path = file("versions.yml")

    script:
    // Derived from task.memory, not hardcoded: leaves headroom below the task limit for
    // the JVM's non-heap overhead, and lets the heap actually grow on a retry now that
    // `memory` scales with task.attempt.
    def xmx = "${(task.memory.giga * 0.8) as int}g"
    """
    set -euo pipefail
    /usr/bin/java -Xmx${xmx} \\
        -Djava.io.tmpdir="\$PWD" \\
        -jar /opt/local/squirls/squirls-cli-1.0.0.jar \\
        annotate-vcf \\
        -t ${task.cpus} \\
        --output-format vcf \\
        "${config}" \\
        "${jannovar}" \\
        "${vcf.name}" \\
        "${tag}.spliceai.squirls"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    squirls: \$(java -jar /opt/local/squirls/squirls-cli-1.0.0.jar --version 2>&1 | head -n1 | sed 's/^squirls v//')
    java: \$(java -version 2>&1 | head -n1 | sed 's/.*"\\(.*\\)".*/\\1/')
END_VERSIONS
    """
}
process compress_squirls {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 10.GB * task.attempt }

    input:
    interval: String
    tag: String
    vcf: Path

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        vcf: file("${vcf.name}.gz"),
        tbi: file("${vcf.name}.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    bgzip "${vcf.name}"
    tabix -p vcf "${vcf.name}.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^tabix (htslib) //')
END_VERSIONS
    """
}
workflow squirls {
    take:
    vcfs: Channel<Record>
    config: Path
    jannovar: Path

    main:
    // multiMap: one traversal, fields stay visibly one row -- see prepare_svs.nf.
    in_parts = vcfs.multiMap { r ->
        interval: r.interval
        tag: r.tag
        vcf: r.vcf
        tbi: r.tbi
    }
    squirls_res = run_squirls(in_parts.interval, in_parts.tag, in_parts.vcf, in_parts.tbi,
                              config, jannovar)

    sq_parts = squirls_res.out.multiMap { r ->
        interval: r.interval
        tag: r.tag
        vcf: r.vcf
    }
    compress_res = compress_squirls(sq_parts.interval, sq_parts.tag, sq_parts.vcf)

    emit:
    compressed: Channel<Record> = compress_res.out
    versions: Channel<Path> = squirls_res.versions.mix(compress_res.versions)
}
