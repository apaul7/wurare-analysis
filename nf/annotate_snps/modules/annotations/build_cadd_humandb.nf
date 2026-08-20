nextflow.enable.types = true

// Builds the ANNOVAR-format CADD database once per run, from the merged genome-wide CADD
// table emitted by merge_cadd.
//
// Extracted out of wf's per-interval run_annovar (annovar.nf:56-82), where it sat inside
// the interval fan-out: `cadd_tsv` is a singleton, identical for every interval, so the
// reformat and index_annovar.pl ran ~24x on a whole-genome run producing byte-identical
// output each time. index_annovar.pl over a genome-wide CADD table is not cheap.
// Same run-once/broadcast shape as upstream's build_vcfanno_conf.
//
// The awk coordinate reformat is carried over verbatim. What it does:
//
// CADD reports the score in column 6 and uses VCF's anchored-base convention for indels.
// ANNOVAR's generic filter format wants the score in column 5, an explicit start/end pair,
// and "-" on whichever side of an indel is empty. So for an insertion
//
//     chr2  8778635  T  TACGAAGAGGG  0.817388  9.566
//
// the shared leading base T is dropped, ref becomes "-", and the row becomes
//
//     chr2  8778635  8778635  -  ACGAAGAGGG  9.566  9.566
//
// Deletions additionally shift start by +1, since the anchor base is not part of the
// deleted sequence. Left alone deliberately -- it is correct and in production, and
// rewriting working indel coordinate arithmetic is not what this port is for.
process build_cadd_humandb {
    container 'apaul7/analysis:1.2.0'
    cpus 1
    memory { 20.GB * task.attempt }

    input:
    tuple(cadd_tsv: Path, cadd_tbi: Path)
    annovar_dir: Path

    output:
    out: Tuple<Path, Path> = tuple(
        file("humandb/hg38_CADDv1.6.txt"),
        file("humandb/hg38_CADDv1.6.txt.idx")
    )
    versions: Path = file("versions.yml")

    script:
    """
    set -euo pipefail
    mkdir humandb

# need to move cadd_phred to col 5
# format:
# chr start stop ref alt score other
# for indels need to use "-" as ref or alt.
zgrep -v "^#" "${cadd_tsv.name}" | \\
    awk -F"\\t" '{OFS="\\t";
\$5=\$6;
if(length(\$3)<length(\$4) && substr(\$4,1,1)==substr(\$3,1,1)){
  \$3="-";
  \$4=substr(\$4,2)};
if(length(\$3)>length(\$4) && substr(\$4,1,1)==substr(\$3,1,1)){
  \$2=\$2+1;
  \$3=substr(\$3,2);
  \$4="-"};
\$2=\$2"\\t"\$2-1+length(\$3);
print \$0}' > hg38_CADDv1.6.txt

    "${annovar_dir}"/index_annovar.pl --filetype A --outfile humandb/hg38_CADDv1.6.txt hg38_CADDv1.6.txt

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    annovar: \$( { "${annovar_dir}"/index_annovar.pl 2>&1 || true; } | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -n1 || true )
    CADD: \$( { zcat "${cadd_tsv.name}" || true; } | head -n1 | sed 's/^##CADD [^-]*-\\(v[0-9.]*\\).*/\\1/' )
END_VERSIONS
    """
}
