// Copied from an upstream pipeline and locally modified; upstream is not tracked.
//
// NO LONGER VERBATIM (beyond the container tag pin, which is repo-wide policy rather than a
// divergence). Local to this repo, and NOT an upstream defect:
//   1. merge_cadd takes `chr_prefix` and re-adds the contig prefix conditionally. Upstream
//      re-adds "chr" unconditionally, which is right for a callset that is always GATK-named;
//      the conditional only matters for the Ensembl-named inputs this pipeline tolerates,
//      where an unconditional prefix yields a table ANNOVAR cannot join -- blank column, no
//      error. The chrM <-> MT special case in both directions came FROM upstream
//      and is composed with the prefix check here.
// Upstream's defects too:
//   2. run_cadd's `input:` takes scalars, not a record -- see nf/shared/modules/other/tools.nf.
//   3. toSortedList, not toList, before merge_cadd -- see the workflow body.
//   4. Versions are read from the tool rather than hardcoded to a literal.
// Local additions upstream has no counterpart for -- the prescored top-up:
//   5. extract_prescored + the workflow's `prescored` take: a prescored CADD table's rows are
//      copied per interval and only the remainder goes to CADD.sh.
//   6. run_cadd guards against a variant-free input, which top-up makes possible (an interval
//      the table fully covers) and CADD.sh errors out on.
//   7. merge_cadd sorts the merged body: prescored fragments and run_cadd outputs interleave
//      within a chromosome, and tabix refuses an unsorted table.
// Also: no `stub:` blocks (this repo tests by running the real DAG against fake tools).
nextflow.enable.types = true
process run_cadd {
    container 'apaul7/docker-cadd:v1.6'
    cpus 10
    memory { 50.GB * task.attempt }

    // Scalars, not a record -- see nf/shared/modules/other/tools.nf.
    input:
    tag: String
    vcf: Path
    tbi: Path
    data_dir: Path

    output:
    out:      Path = file("${tag}.caddv1.6.out.tsv.gz")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # Strip the `chr` prefix CADD's GRCh38 data does not use, and rename the mitochondrion.
    #
    # Stripping `chr` alone is not enough: CADD v1.6 GRCh38 is Ensembl-named, where the
    # mitochondrial contig is MT, not M. CADD's own code agrees -- see
    # src/scripts/lib/AnalysisLib.py, which normalises `M` to `MT` before looking a contig up
    # in Ensembl's files. A bare prefix strip sends `M`, which its MT-keyed lookups miss, and
    # the miss is silent because a CADD fragment with no rows for an interval is
    # indistinguishable from an interval that legitimately has none (see merge_cadd's
    # `|| true`).
    #
    # awk on field 1, not sed on the line: sed's `\\t` is a GNU extension and this image ships
    # mawk, and an unanchored `^chrM` would also rewrite any contig merely starting with those
    # characters. The inverse mapping lives in merge_cadd and must stay in step with this one
    # -- send MT and read MT back, or the pair silently stops matching.
    zgrep "^#" "${vcf.name}" > nochr.vcf
    # `|| true`: a variant-free VCF is a legitimate input when a prescored table covered the
    # whole interval (extract_prescored), and zgrep exits 1 on no-match, which pipefail would
    # turn into an abort.
    { zgrep -v "^#" "${vcf.name}" || true; } \\
      | awk 'BEGIN{FS=OFS="\\t"} {sub(/^chr/, "", \$1); if (\$1 == "M") \$1 = "MT"; print}' \\
      >> nochr.vcf

    # CADD.sh errors out on a VCF with no variants (the same fact that makes list_chroms
    # derive intervals from the input rather than a static chr1..chrY list). An interval the
    # prescored table fully covered still owes merge_cadd an output, so emit a table with a
    # header and no rows -- merge_cadd reads the merged header from its first file, and
    # build_cadd_humandb parses the version out of the ##CADD line, so both lines matter.
    if grep -qv "^#" nochr.vcf; then
        "${data_dir}"/CADD.sh \\
            -g GRCh38 \\
            -v v1.6 \\
            -c ${task.cpus} \\
            -o "${tag}.caddv1.6.out.tsv.gz" \\
            nochr.vcf
    else
        printf '##CADD GRCh38-v1.6 (no unscored variants in this interval)\\n#Chrom\\tPos\\tRef\\tAlt\\tRawScore\\tPHRED\\n' \\
          | gzip -c > "${tag}.caddv1.6.out.tsv.gz"
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    CADD: \$( { "${data_dir}"/CADD.sh -h 2>&1 || true; } | head -n1 | sed 's/.*-- CADD version /v/' )
END_VERSIONS
    """
}

process merge_cadd {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 5.GB * task.attempt }

    input:
    tsvs: List<Path>
    // "chr" or "". run_cadd strips a leading "chr" conditionally; this re-adds it
    // symmetrically. Derived from the intervals -- see the file header.
    chr_prefix: String

    output:
    out:      Tuple<Path, Path> = tuple(file('caddv1.6.out.tsv.gz'), file('caddv1.6.out.tsv.gz.tbi'))
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    zgrep "^#" ${tsvs.head().name} > caddv1.6.out.tsv
    for t in ${tsvs.join(" ")}; do
        # `|| true`: a per-interval CADD fragment with a header but no data rows is a
        # legitimate result, and zgrep exits 1 on no-match. Under `set -o pipefail` that
        # would abort the whole merge over an empty interval.
        #
        # The inverse of run_cadd's mapping, and it has to special-case MT: prefixing
        # everything would turn CADD's MT into chrMT, which matches no contig in the input
        # VCF, so ANNOVAR's join would miss every mitochondrial variant just as surely as
        # dropping them did. Keep this in step with run_cadd above.
        #
        # Nothing is rewritten when the callset is Ensembl-named (chr_prefix ""), because
        # run_cadd's strip was then a no-op and MT is already the input's own name.
        { zgrep -v "^#" \$t || true; } \\
          | awk -v p='${chr_prefix}' 'BEGIN{FS=OFS="\\t"} {if (p != "") {if (\$1 == "MT") \$1 = "chrM"; else \$1 = p \$1} print}'
    done \\
      | LC_ALL=C sort -k1,1 -k2,2n -s >> caddv1.6.out.tsv
    # The sort exists for the prescored top-up: extract_prescored's fragment and run_cadd's
    # output cover the SAME chromosome, so concatenating them leaves positions unsorted
    # within a contig and `tabix` below refuses the file. Without a prescored table it is a
    # near-no-op -- the filename-sorted input list already arrives in lexical contig order
    # with each file internally position-sorted -- and `-s` keeps it deterministic either way.
    bgzip caddv1.6.out.tsv
    tabix -s 1 -b 2 -e 2 caddv1.6.out.tsv.gz

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^tabix (htslib) //')
END_VERSIONS
    """
}

// Splits an interval VCF against a prescored CADD table: matches become a per-interval CADD
// output fragment copied from the table, everything else goes on to run_cadd. The CADD
// counterpart of spliceai.nf's add_precomputed, except the scored half is a TSV, not a VCF.
process extract_prescored {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 10.GB * task.attempt }

    input:
    interval: String
    tag: String
    vcf: Path
    tbi: Path
    prescored: Path
    // A real input so it is staged next to `prescored` and pysam can seek.
    prescored_tbi: Path

    output:
    out: Record = record(
        tag: tag,
        unscored_vcf: file("${tag}.unscored.vcf.gz"),
        unscored_tbi: file("${tag}.unscored.vcf.gz.tbi"),
        // Named so it sorts next to run_cadd's ${tag}.caddv1.6.out.tsv.gz in merge_cadd's
        // filename-sorted input list.
        fragment: file("${tag}.prescored.caddv1.6.out.tsv.gz")
    )
    versions: Path = file("versions.yml")

    // add_cadd_scores.py lives in ../../bin/, which Nextflow puts on PATH automatically.
    script:
    """
    set -euo pipefail
    add_cadd_scores.py \\
        --input "${vcf.name}" \\
        --prescored "${prescored}" \\
        --region "${interval}" \\
        --scored-tsv "${tag}.prescored.caddv1.6.out.tsv" \\
        --unscored "${tag}.unscored.vcf.gz"
    tabix -p vcf "${tag}.unscored.vcf.gz"
    bgzip "${tag}.prescored.caddv1.6.out.tsv"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/^Python //')
    pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^tabix (htslib) //')
END_VERSIONS
    """
}

workflow cadd {
    take:
    vcfs: Channel<Record>
    data_dir: Path
    // Empty list, or [tsv_gz, tbi]: a prescored CADD table whose scores are copied so only
    // the remainder is scored. Same absence convention as the subworkflow's precomputed_cadd.
    prescored: List<Path>

    main:
    // The intervals are list_chroms' output, i.e. the input VCF's own CHROM values, so they
    // are the authority on whether this callset uses a "chr" prefix. Taken from the first
    // interval: a VCF mixing the two conventions is malformed and bcftools would have
    // rejected it long before here. Derived BEFORE the multiMap below, which consumes the
    // record channel.
    chr_prefix = vcfs.first().map { r -> r.interval.startsWith('chr') ? 'chr' : '' }

    // multiMap: one traversal, and the fields stay visibly one row. See prepare_svs.nf --
    // separate .map reads would also work, a queue channel is broadcast to every consumer.
    if (prescored) {
        in_parts = vcfs.multiMap { r ->
            interval: r.interval
            tag: r.tag
            vcf: r.vcf
            tbi: r.tbi
        }
        extracted = extract_prescored(in_parts.interval, in_parts.tag,
                                      in_parts.vcf, in_parts.tbi,
                                      prescored[0], prescored[1])
        ex_parts = extracted.out.multiMap { r ->
            tag: r.tag
            vcf: r.unscored_vcf
            tbi: r.unscored_tbi
        }
        cadd_result = run_cadd(ex_parts.tag, ex_parts.vcf, ex_parts.tbi, data_dir)
        // The fragments join run_cadd's outputs in ONE merge: both are per-interval CADD
        // tables in native naming, and merge_cadd's body sort interleaves them.
        tsvs = cadd_result.out.mix(extracted.out.map { r -> r.fragment })
        step_versions = cadd_result.versions.mix(extracted.versions)
    } else {
        parts = vcfs.multiMap { r ->
            tag: r.tag
            vcf: r.vcf
            tbi: r.tbi
        }
        cadd_result = run_cadd(parts.tag, parts.vcf, parts.tbi, data_dir)
        tsvs = cadd_result.out
        step_versions = cadd_result.versions
    }
    // toSortedList, not toList: collected items arrive in completion order, so both the
    // merged table's row order and this input's hash varied run to run.
    merged = merge_cadd(tsvs.toSortedList { a, b -> a.name <=> b.name },
                        chr_prefix)

    emit:
    cadd_merged: Tuple<Path, Path> = merged.out
    versions: Channel<Path> = step_versions.mix(merged.versions)
}
