nextflow.enable.types = true
include { Reference } from "../../../shared/modules/types.nf"

// Rewritten from wf/modules/annotations/annovar.nf, which ran all six steps as one
// monolithic per-interval process with numbered echo progress markers. The ANNOVAR
// commands and their flags carry over character-for-character; only the structure around
// them changed. Each step is now its own process, so the Nextflow trace *is* the progress
// log and the echoes are gone -- and each step gets its own memory profile and its own
// retry instead of sharing one flat 50 GB.

// ANNOVAR ships no --version flag: each perl script prints a `Version: $Date: ...` banner
// on bare invocation and exits non-zero, hence the `|| true`. Probed rather than pinned --
// an ANNOVAR install and its humandb drift over time, and a hardcoded version string would
// go stale without anyone noticing.
def annovar_version(annovar_dir, tool) {
    return "\$( { \"${annovar_dir}\"/${tool} 2>&1 || true; } | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -n1 || true )"
}

// ANNOVAR requires every database to live in one directory, but the CADD db is built at
// runtime by build_cadd_humandb while the rest ship with the install -- so both are
// symlinked into a task-local ./humandb/ instead of writing into the shared, read-only
// install directory.
//
// ${PWD} is load-bearing. `ln -s` resolves a relative target against the *link's own*
// directory (./humandb/), not the task working directory, so a relative target here
// produces links that silently dangle.
//
// One definition, called by both reduce_variants and table_annovar -- the split into two
// processes duplicates the farm at runtime, not in source.
def stage_humandb(annovar_dir, cadd_txt, cadd_idx) {
    return """mkdir humandb
    ln -s "\${PWD}/${annovar_dir}/humandb/"*.txt ./humandb/
    ln -s "\${PWD}/${annovar_dir}/humandb/"*.txt.idx ./humandb/
    ln -s "\${PWD}/${annovar_dir}/humandb/"*.fa ./humandb/
    ln -s "\${PWD}/${cadd_txt.name}" ./humandb/
    ln -s "\${PWD}/${cadd_idx.name}" ./humandb/"""
}

process normalize_vcf {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 10.GB * task.attempt }

    input:
    interval: String
    tag: String
    vcf: Path
    tbi: Path
    ref: Reference

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        vcf: file("${tag}.av.vcf")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # Piped rather than writing a tmp.vcf and rm-ing it; -Ou keeps the intermediate as
    # uncompressed BCF on the pipe.
    bcftools norm -m-both -Ou "${vcf.name}" \\
        | bcftools norm --check-ref s --fasta-ref "${ref.fa}" -o "${tag}.av.vcf"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// Steps 02 and 03 of the original stay fused deliberately: the per-sample loop reads the
// `-outfile sample -allsample` scratch files convert2annovar.pl drops in the same working
// directory. Splitting them would mean declaring an unknown-at-compile-time set of
// intermediate files as outputs.
process make_avinput {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 10.GB * task.attempt }

    input:
    interval: String
    tag: String
    vcf: Path
    annovar_dir: Path

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        avinput: file("${tag}.sorted.avinput")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    "${annovar_dir}"/convert2annovar.pl \\
        -format vcf4 \\
        "${vcf.name}" \\
        -outfile sample \\
        --includeinfo \\
        -allsample

    # Globbed rather than recovered via `ls | tr | grep sample | sed | sed`, which also
    # matched any sample whose own name contained "sample".
    for f in sample.*.avinput; do
        s="\${f#sample.}"
        s="\${s%.avinput}"
        awk -v s="\$s" '{print \$0"\\t"s}' "\$f" > "\$f.tag"
    done
    cat ./sample.*.avinput.tag > all_sample.tag.avinput
    bedtools sort -chrThenSizeA -i all_sample.tag.avinput > "${tag}.sorted.avinput"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    annovar: ${annovar_version(annovar_dir, 'convert2annovar.pl')}
    bedtools: \$(bedtools --version 2>&1 | head -n1 | sed 's/^bedtools //')
END_VERSIONS
    """
}

// Lifts the SpliceAI and SQUIRLS INFO fields carried through from the upstream VCF into
// their own avinput columns. Its own process because it is a distinct tool (python) with a
// distinct purpose from the ANNOVAR perl steps either side of it.
process add_splice_scores {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 10.GB * task.attempt }

    input:
    interval: String
    tag: String
    avinput: Path
    splice_scores_script: Path

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        avinput: file("${tag}.splice_scores.avinput")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    python3 "${splice_scores_script}" \\
        --input "${avinput.name}" \\
        --output "${tag}.splice_scores.avinput"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/^Python //')
END_VERSIONS
    """
}

// variants_reduction.pl runs one numbered step per --protocol entry and writes
// <outfile>.stepN.* as it goes. Only .step2.varlist is declared as output, and that "2" is
// tied to --protocol/--operation being exactly two entries
// (gnomad40_genome,gnomad40_exome / f,f). Change the reduction protocol list and this
// filename changes with it -- otherwise the handoff to table_annovar breaks silently.
process reduce_variants {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 50.GB * task.attempt }

    input:
    interval: String
    tag: String
    avinput: Path
    annovar_dir: Path
    tuple(cadd_txt: Path, cadd_idx: Path)

    output:
    out: Record = record(
        interval: interval,
        tag: tag,
        varlist: file("${tag}.reduced.avinput.step2.varlist")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    ${stage_humandb(annovar_dir, cadd_txt, cadd_idx)}

    "${annovar_dir}"/variants_reduction.pl \\
        "${avinput.name}" ./humandb \\
        --buildver hg38 \\
        --protocol gnomad40_genome,gnomad40_exome \\
        --operation f,f \\
        --aaf_threshold 0.01 \\
        --remove \\
        --outfile "${tag}.reduced.avinput"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    annovar: ${annovar_version(annovar_dir, 'variants_reduction.pl')}
END_VERSIONS
    """
}

process table_annovar {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 50.GB * task.attempt }

    input:
    interval: String
    tag: String
    varlist: Path
    annovar_dir: Path
    tuple(cadd_txt: Path, cadd_idx: Path)
    // Resolved in main.nf (CLINVAR already substituted) so the humandb preflight and this
    // command line are guaranteed to be the same list. Positions must correspond.
    protocols: String
    operations: String
    omim_xref: Path
    out_prefix: String

    output:
    out:      Path = file("${out_prefix}.${tag}.hg38_multianno.txt")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    ${stage_humandb(annovar_dir, cadd_txt, cadd_idx)}

    "${annovar_dir}"/table_annovar.pl \\
        "${varlist.name}" ./humandb/ \\
        --buildVer hg38 \\
        --out "${out_prefix}.${tag}" \\
        --remove \\
        --protocol ${protocols} \\
        --operation ${operations} \\
        --nastring . \\
        --otherinfo \\
        --polish \\
        -xref "${omim_xref}"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    annovar: ${annovar_version(annovar_dir, 'table_annovar.pl')}
END_VERSIONS
    """
}

// Text-level merge of the per-interval multianno fragments: header once, then tail-concat.
// Same shape as upstream's merge_annotations -- the per-interval files are never
// recombined as VCFs.
process merge_annovar {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 5.GB * task.attempt }

    input:
    tsvs: List<Path>
    out_prefix: String

    output:
    out:      Path = file("${out_prefix}.hg38_multianno.tsv")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    head -n1 "${tsvs.head().name}" > "${out_prefix}.hg38_multianno.tsv"
    for t in ${tsvs.join(" ")}; do
        tail -n+2 \$t >> "${out_prefix}.hg38_multianno.tsv"
    done

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    coreutils: \$(head --version 2>&1 | head -n1 | grep -oE '[0-9]+\\.[0-9]+' | head -n1 || true)
END_VERSIONS
    """
}

workflow annovar {
    take:
    vcfs: Channel<Record>
    cadd_humandb: Tuple<Path, Path>
    cohort: String
    data_type: String
    ref: Reference
    clinvar_date: String
    // Already CLINVAR-substituted and preflighted in main.nf.
    protocols: String
    operations: String
    omim_xref: Path
    annovar_dir: Path
    splice_scores_script: Path

    main:
    // Computed once, here, and passed down as a String. wf derived the final filename from
    // the first list element's basename truncated at its first dot, making it an emergent
    // property of an arbitrary list element that only held while cohort and data_type
    // themselves contained no dots.
    def run_date = clinvar_date.tokenize('_')[1]
    def out_prefix = "${cohort}_${data_type}_${run_date}".toString()

    // Records destructured before each process sees them, so -resume can hash by value --
    // see nf/shared/modules/other/tools.nf. They stay records in every `output:`.
    vcf_parts = vcfs.multiMap { r ->
        interval: r.interval
        tag: r.tag
        vcf: r.vcf
        tbi: r.tbi
    }
    normalized = normalize_vcf(vcf_parts.interval, vcf_parts.tag,
                               vcf_parts.vcf, vcf_parts.tbi, ref)

    norm_parts = normalized.out.multiMap { r ->
        interval: r.interval
        tag: r.tag
        vcf: r.vcf
    }
    avinput = make_avinput(norm_parts.interval, norm_parts.tag, norm_parts.vcf, annovar_dir)

    av_parts = avinput.out.multiMap { r ->
        interval: r.interval
        tag: r.tag
        avinput: r.avinput
    }
    scored = add_splice_scores(av_parts.interval, av_parts.tag, av_parts.avinput,
                               splice_scores_script)

    scored_parts = scored.out.multiMap { r ->
        interval: r.interval
        tag: r.tag
        avinput: r.avinput
    }
    reduced = reduce_variants(scored_parts.interval, scored_parts.tag, scored_parts.avinput,
                              annovar_dir, cadd_humandb)

    red_parts = reduced.out.multiMap { r ->
        interval: r.interval
        tag: r.tag
        varlist: r.varlist
    }
    tabled = table_annovar(red_parts.interval, red_parts.tag, red_parts.varlist,
                           annovar_dir, cadd_humandb, protocols, operations,
                           omim_xref, out_prefix)
    // toSortedList, not toList: collected items arrive in completion order, so both the
    // merged TSV's interval order and this list's hash varied run to run.
    merged     = merge_annovar(tabled.out.toSortedList { a, b -> a.name <=> b.name },
                               out_prefix)

    emit:
    annovar_tsv: Path = merged.out
    versions: Channel<Path> = normalized.versions
        .mix(avinput.versions, scored.versions, reduced.versions, tabled.versions, merged.versions)
}
