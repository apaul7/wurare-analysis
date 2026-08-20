nextflow.enable.types = true
include { MergedVcf } from '../merge/svdb_merge.nf'

// Phase 6: the Talos tail.
//
// Talos accepts one schema and nothing else, which is why this is a tail hanging off the
// shared backbone rather than the pipeline's shape. Everything before it serves every
// consumer; only this stage is Talos-specific, and it is optional.
//
// Deliberate simplification against the Talos plan's stage 7. That stage derives external AF
// through GATK-SV's own toolchain -- bedtools closest, two R scripts, and
// annotate_external_af_modify_vcf.py, all extracted from a pinned sv-pipeline image. This
// pipeline already has gnomAD frequency on the VCF from Phase 4's svdb query, so the tail
// renames that instead of pulling in a second heavyweight image to recompute the same
// quantity by a different route. Talos needs the FIELD, not the derivation.
//
// What that costs, stated rather than buried: `{GNOMAD_POP}_sv_SVID` has no value. svdb
// query returns an occurrence count and a frequency, not the identifier of the matching
// gnomAD variant. Talos uses SVID for reporting rather than filtering, so a null there is a
// degraded report and not a wrong one -- if SVID turns out to matter, stage 7's bedtools/R
// route is what supplies it.
//
// The field is still DECLARED in the header, and that distinction is the whole point.
// Talos's rearrange_annotations() tolerates a missing ALGORITHMS, STATUS, CHR2 and END2 and
// nothing else; the rest, SVID included, are `mt.info[...]` struct accesses, and hail builds
// that struct from the VCF header rather than from the records. Omitting the field is a hail
// error before any filtering runs -- Talos does not start. Declaring it empty is the degraded
// report. assets/talos_schema.awk backfills every such header.

// PREDICTED_LOF. Talos hard gates on this: a variant without at least one entry is dropped
// entirely, so this process is not optional within the tail even though the tail is.
process svannotate {
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_0'
    cpus 2
    memory { 16.GB * task.attempt }

    input:
    label: String
    vcf: Path
    tbi: Path
    protein_coding_gtf: Path
    noncoding_bed: Path
    // Rewrites sequence-resolved ALTs to symbolic ones. Not cosmetic: SVAnnotate throws
    // IllegalArgumentException on a literal ALT and takes the whole tail down with it, and
    // Manta emits literal ALTs as a matter of course.
    symbolic_awk: Path

    output:
    out: Path = file("${label}.svannotate.vcf")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # gzip rather than bcftools: this is the gatk4 container, which carries no bcftools, and
    # adding a second image to decompress one file buys nothing.
    gzip -dc "${vcf.name}" | awk -f "${symbolic_awk.name}" > svannotate_input.vcf

    gatk SVAnnotate \\
        -V svannotate_input.vcf \\
        --protein-coding-gtf "${protein_coding_gtf.name}" \\
        --non-coding-bed "${noncoding_bed.name}" \\
        -O "${label}.svannotate.vcf"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    gatk4: \$(gatk --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1)
END_VERSIONS
    """
}

// Sex-stratified AF, from the PED. Talos filters the recessive path on AF_MALE/AF_FEMALE.
// bcftools +fill-tags with a sample-group file gives these directly, which avoids depending
// on GATK-SV's compute_AFs.py and the image it lives in.
process sex_stratified_af {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 8.GB * task.attempt }

    input:
    label: String
    vcf: Path
    ped: Path

    output:
    out: Path = file("${label}.sexaf.vcf")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # PED column 5 is sex: 1 male, 2 female, 0 unknown. Unknown samples are left out of both
    # groups rather than guessed into one -- a wrong sex assignment shifts an AF that a
    # recessive filter then acts on.
    awk 'BEGIN{FS="[ \\t]+"} \$0 !~ /^#/ && NF >= 5 {
             if (\$5 == 1) print \$2, "MALE";
             else if (\$5 == 2) print \$2, "FEMALE"
         }' OFS="\\t" "${ped.name}" > groups.txt

    # fill-tags declares only the groups the PED actually contains, so a single-sex cohort
    # gets AF_MALE and no AF_FEMALE -- and Talos, seeing AF_MALE, reads AF_FEMALE off the
    # same struct and raises. talos_schema downstream backfills whichever header is absent;
    # nothing here fabricates the value.
    if [ -s groups.txt ]; then
        bcftools +fill-tags "${vcf.name}" -Ov -o "${label}.sexaf.vcf" -- -S groups.txt -t AF
    else
        echo "WARNING: no sexed samples in the PED; AF_MALE/AF_FEMALE not computed" >&2
        cp "${vcf.name}" "${label}.sexaf.vcf"
    fi

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}

// Rename and derive the rest of the schema, then report what is present. The report is
// published: Talos drops what it cannot read without saying so, and this is the only place
// that says so.
process talos_schema {
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    label: String
    vcf: Path
    schema_awk: Path
    keep_awk: Path
    check_awk: Path
    gnomad_pop: String
    joint: Boolean

    output:
    out: MergedVcf = new MergedVcf(
        label: label,
        joint: joint,
        vcf: file("${label}.talos.vcf.gz"),
        tbi: file("${label}.talos.vcf.gz.tbi")
    )
    report: Path = file("${label}.talos_fields.txt")
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    awk -v GNOMAD_POP="${gnomad_pop}" -f "${schema_awk.name}" "${vcf.name}" > schema.vcf

    # Strip INFO down to Talos's read set -- PRPOS/PREND and the rest of the input
    # callers' INFO baggage otherwise ride into the handoff file. After the schema awk,
    # so the derived fields (SOFT_FILTERS, N_HET, ...) exist to be kept.
    # talos_keep_info.awk exits nonzero on an empty keep-list, which fails this
    # assignment under set -e rather than handing bcftools a bare "^" that would strip
    # every INFO field.
    KEEP=\$(awk -v GNOMAD_POP="${gnomad_pop}" -f "${keep_awk.name}" schema.vcf)
    bcftools annotate -x "\$KEEP" schema.vcf -Oz -o "${label}.talos.vcf.gz"
    bcftools index -t "${label}.talos.vcf.gz"

    bcftools view "${label}.talos.vcf.gz" \\
        | awk -v GNOMAD_POP="${gnomad_pop}" -f "${check_awk.name}" \\
        > "${label}.talos_fields.txt"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}
