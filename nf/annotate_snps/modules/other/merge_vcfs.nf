nextflow.enable.types = true

// Recombines the per-interval annotated VCFs into one cohort VCF.
//
// The upstream pipeline's annotate_variants deliberately never does this -- there the per-interval
// VCFs are terminal and only the text-level annotations.tsv is merged. Here the merged VCF
// is a published artifact in its own right: it carries the SpliceAI and SQUIRLS INFO fields,
// which are the expensive part of a run, so re-running ANNOVAR against updated databases
// needs only this file plus the merged CADD table -- not another CADD/SpliceAI/SQUIRLS pass.
//
// NOT a fully annotated VCF. CADD scores and every ANNOVAR field are absent: CADD is emitted
// as its own TSV and joined by ANNOVAR at the humandb level, and ANNOVAR writes TSV, never
// VCF. Unlike upstream there is no vcfanno step to fold CADD into INFO. The filename
// says spliceai_squirls rather than "annotated" for exactly that reason.
//
// No run-date component in the name: the contents depend on SpliceAI/SQUIRLS only, not on
// clinvar_date, so borrowing the multianno TSV's date would imply a coupling that isn't real.
//
// --allow-overlaps is used with indexed inputs so bcftools orders records by the index rather
// than trusting the order the interval channel happened to emit. The intervals are distinct
// contigs and cannot actually overlap today, but relying on channel ordering for genomic sort
// order would be a silent correctness bug the first time that stops holding.
process merge_vcfs {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 10.GB * task.attempt }

    input:
    vcfs: List<Path>
    tbis: List<Path>
    out_prefix: String

    output:
    out: Tuple<Path, Path> = tuple(
        file("${out_prefix}.spliceai_squirls.vcf.gz"),
        file("${out_prefix}.spliceai_squirls.vcf.gz.tbi")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    # tbis is never named in the command line, on purpose: declaring it as an input is what
    # stages each .tbi beside its .vcf.gz so --allow-overlaps can read the indexes.
    bcftools concat --allow-overlaps -Oz -o "${out_prefix}.spliceai_squirls.vcf.gz" ${vcfs.join(" ")}
    bcftools index -t "${out_prefix}.spliceai_squirls.vcf.gz"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
END_VERSIONS
    """
}
