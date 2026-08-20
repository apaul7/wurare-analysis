// Copied from the upstream pipeline this was ported from.
//
// NO LONGER VERBATIM. Two local hunks: the container is pinned to a tag rather than left at
// :latest, and the `stub:` block is absent (this repo tests by running the real DAG against
// fake tools instead -- see tests/run.sh). Everything else is unchanged.
nextflow.enable.types = true

// Collects the per-task versions.yml fragments every process emits into one file published
// at <outdir>/pipeline_info/software_versions.yml.
//
// Each fragment is a YAML block keyed by process name, e.g.
//
//   prep_and_merge_sample_vcf:
//     bcftools: 1.21
//     bedtools: 2.31.1
//
// Fragments are deduplicated: a process that ran 200 per-interval tasks reports its tool
// versions once, not 200 times. Sorted so the published file is byte-identical across runs
// with the same software.
//
// Call this on a channel of fragment Paths via mergeVersions() below -- not by calling
// dump_versions directly -- since the merge step matters as much as the process itself.
//
// This is two steps, not one process, because of two Nextflow behaviors confirmed by
// actually running a pipeline (not assumed from docs):
//
//   1. `exec:` blocks do not stage process inputs. A Path passed into one is a
//      `nextflow.processor.TaskPath`, which reports the correct FileSystemProvider but
//      isn't the concrete platform Path class the JDK's native filesystem implementation
//      requires -- any real I/O on it (`.text`, `.copyTo()`, even `.toAbsolutePath()`)
//      throws (ProviderMismatchException / UnsupportedOperationException). This is what
//      broke the original single-process exec: version on a real compute2 run. Reading
//      fragment *content* has to happen in a `script:` block, where Nextflow does stage
//      inputs as real files.
//   2. Every process emits its fragment as the literal filename "versions.yml". Handing
//      that many identically-named files to one process's `List<Path>` input does not
//      auto-disambiguate on staging -- confirmed: with 3 same-named inputs, only 1 symlink
//      is created and the other 2 are silently dropped, silently losing data. Typed
//      syntax also doesn't accept the classic DSL2 `stageAs:` input modifier that would
//      normally solve this.
//
// The fix: merge with `collectFile()` first -- a core Nextflow channel operator that reads
// file content directly, outside any process's staging, so neither problem applies -- then
// hand the one resulting file (now a single, uniquely-named input) to a `script:` process
// that dedupes/sorts the blank-line-separated blocks.
workflow mergeVersions {
    take:
    versions_ch: Channel<Path>

    main:
    raw = versions_ch.collectFile(name: 'raw_versions.yml', newLine: true)
    merged = dump_versions(raw)

    emit:
    merged
}

process dump_versions {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 1.GB * task.attempt }

    input:
    raw: Path

    output:
    file("software_versions.yml")

    script:
    """
    set -euo pipefail
    cat > dedupe_versions.py <<'PYEOF'
#!/usr/bin/env python3
import sys
with open(sys.argv[1]) as fh:
    blocks = [b.strip() for b in fh.read().split("\\n\\n")]
blocks = sorted(set(b for b in blocks if b))
with open("software_versions.yml", "w") as out:
    out.write("\\n\\n".join(blocks) + "\\n")
PYEOF
    python3 dedupe_versions.py "${raw.name}"
    """
}
