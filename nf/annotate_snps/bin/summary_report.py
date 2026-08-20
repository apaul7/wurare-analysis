#!/usr/bin/env python3
"""One tidy TSV summarizing an annotate_snps run: section <tab> metric <tab> sample <tab> value.

    summary_report.py --stats bcftools_stats.txt --multianno X.hg38_multianno.tsv \
        --cohort C --data-type T --run-date D --clinvar-date CD \
        [--samples-tsv F] [--ploidy F] [--ancestry F] [--pairs F] --out summary.tsv

Cohort-level rows use sample ".". Sections, in fixed output order:

  run          cohort/data_type/run_date/clinvar_date/n_samples
  cohort       whole-VCF counts from `bcftools stats` (SN + TSTV)
  sample       per-sample counts from the PSC block
  rare_subset  the multianno TSV -- which holds only the RARE subset: reduce_variants
               filters to gnomAD AAF < 0.01 before table_annovar ever runs
  flags        thresholded counts over the same rare subset
  qc           somalier-derived sex/ancestry/relatedness, when those files were given

Tolerant by design: this is a report, and a report must never kill the run that produced
the data it summarizes. A missing file, an empty stats output, or a column absent from a
custom --annovar_protocols list drops the affected rows with a note on stderr and exits 0.
Columns are located by NAME, never by position (the check_somalier.awk discipline) -- but
where that file fails hard on a rename, this one skips, because nothing downstream consumes
these numbers.

Multianno rows are (variant x carrier-sample) pairs -- make_avinput -allsample writes one
line per carrying sample -- so every count here is over UNIQUE variants (first-seen dedup on
Chr:Start:End:Ref:Alt) except the row count that says so in its name.

Output row order is fixed and samples are sorted, so the file is byte-stable across runs
and does not defeat -resume for anything that might one day read it.
"""

import argparse
import re
import sys

# Thresholds for the flags section. Constants, not flags: promote to params only if someone
# actually needs to tune them. (The rare cutoff itself, 0.01, lives upstream in
# reduce_variants.)
CADD_MIN = 20.0
SPLICE_MIN = 0.5

# Same cut points as nf/shared/assets/check_somalier.awk -- see there for why depth ratios,
# not somalier's own sex column, and why the gaps between each pair are deliberately wide.
X_HEMI_MAX = 0.65
X_DIPLOID_MIN = 0.80
Y_ABSENT_MAX = 0.05
Y_PRESENT_MIN = 0.15


def note(msg):
    print(f"summary_report: NOTE {msg}", file=sys.stderr)


def read_lines(path, what):
    """Lines of a TSV, or None (with a note) if absent, unreadable or empty."""
    if not path:
        return None
    try:
        with open(path) as fh:
            lines = fh.read().splitlines()
    except OSError as e:
        note(f"cannot read {what} ({e}); its rows are omitted")
        return None
    if not lines:
        note(f"{what} is empty; its rows are omitted")
        return None
    return lines


def header_index(header, *names):
    """0-based column lookup by name ('#' stripped); first matching candidate, else None."""
    cols = [c.lstrip("#") for c in header.split("\t")]
    for name in names:
        if name in cols:
            return cols.index(name)
    return None


def numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# --- bcftools stats -----------------------------------------------------------------------

# SN key -> metric name. Keys are bcftools' literal text, values ours.
SN_METRICS = [
    ("number of samples:", "n_samples"),
    ("number of records:", "total_records"),
    ("number of SNPs:", "snps"),
    ("number of indels:", "indels"),
    ("number of MNPs:", "mnps"),
    ("number of others:", "others"),
    ("number of multiallelic sites:", "multiallelic_sites"),
    ("number of multiallelic SNP sites:", "multiallelic_snp_sites"),
]

PSC_METRICS = [
    ("nNonRefHom", "n_non_ref_hom"),
    ("nHets", "n_het"),
    (None, "het_hom_ratio"),  # derived below
    ("nSingletons", "n_singletons"),
    ("nMissing", "n_missing"),
    ("nIndels", "n_indels"),
    ("average depth", "mean_depth"),
]


def parse_stats(path):
    """Returns (sn: dict, tstv: str|None, psc: {sample: {colname: value}})."""
    sn, tstv, psc = {}, None, {}
    lines = read_lines(path, "bcftools stats output")
    if lines is None:
        return sn, tstv, psc

    # bcftools names each section's columns in a "# PSC\t[3]sample\t..." header line; map
    # name -> field index rather than trusting positions.
    def section_cols(prefix):
        for line in lines:
            if line.startswith(f"# {prefix}\t") and "[" in line:
                return {re.sub(r"^\[\d+\]", "", f): i
                        for i, f in enumerate(line[2:].split("\t"))}
        return {}

    for line in lines:
        if line.startswith("SN\t"):
            parts = line.split("\t")
            for key, metric in SN_METRICS:
                if len(parts) >= 4 and parts[2] == key:
                    sn[metric] = parts[3]

    tstv_cols = section_cols("TSTV")
    for line in lines:
        if line.startswith("TSTV\t"):
            parts = line.split("\t")
            i = tstv_cols.get("ts/tv")
            if i is not None and i < len(parts):
                tstv = parts[i]

    psc_cols = section_cols("PSC")
    for line in lines:
        if line.startswith("PSC\t"):
            parts = line.split("\t")
            i_sample = psc_cols.get("sample")
            if i_sample is None or i_sample >= len(parts):
                continue
            psc[parts[i_sample]] = {name: parts[i] for name, i in psc_cols.items()
                                    if i < len(parts)}

    if not sn:
        note("bcftools stats output has no SN section; cohort/sample rows are omitted")
    return sn, tstv, psc


# --- multianno ----------------------------------------------------------------------------

def parse_multianno(path, sample_names):
    """One streaming pass. Returns (rare: dict, flags: dict, per_sample: {sample: n})."""
    rare, flags, per_sample = {}, {}, {}
    lines = read_lines(path, "multianno TSV")
    if lines is None:
        return rare, flags, per_sample
    header = lines[0].split("\t")

    key_idx = [header_index(lines[0], n) for n in ("Chr", "Start", "End", "Ref", "Alt")]
    if None in key_idx:
        note("multianno TSV lacks Chr/Start/End/Ref/Alt; rare_subset/flags rows are omitted")
        return rare, flags, per_sample

    i_func = header_index(lines[0], "Func.refGene")
    i_exonic = header_index(lines[0], "ExonicFunc.refGene")
    i_clnsig = header_index(lines[0], "CLNSIG")
    # The runtime-built CADD humandb is headerless, so table_annovar names its column after
    # the protocol; CADD_phred covers a db built with its own header.
    i_cadd = header_index(lines[0], "CADDv1.6", "CADD_phred")
    # Every gnomAD frequency column is named AF (one per protocol; table_annovar does not
    # dedupe). Novel therefore means: absent from every gnomAD set annotated.
    i_afs = [i for i, c in enumerate(header) if c == "AF"]
    # The splice-score columns come from the host-local --annovar_splice_scores_script,
    # whose column names this repo does not control -- matched by substring, or skipped.
    i_spliceai = [i for i, c in enumerate(header) if "spliceai" in c.lower()]
    i_squirls = [i for i, c in enumerate(header) if "squirls" in c.lower()]
    for what, found in [("Func.refGene", i_func is not None),
                        ("ExonicFunc.refGene", i_exonic is not None),
                        ("CLNSIG", i_clnsig is not None),
                        ("CADDv1.6/CADD_phred", i_cadd is not None),
                        ("gnomAD AF", bool(i_afs)),
                        ("spliceai", bool(i_spliceai)),
                        ("squirls", bool(i_squirls))]:
        if not found:
            note(f"multianno TSV has no {what} column; its rows are omitted")

    # The carrier-sample column make_avinput appended lands among the anonymous Otherinfo
    # columns, at a position the host splice-scores script may shift -- so it is recognized
    # by CONTENT: an Otherinfo column whose every value is a VCF sample name. Last match
    # wins (the sample is appended after the includeinfo columns).
    otherinfo = [i for i, c in enumerate(header) if re.fullmatch(r"Otherinfo\d*", c)]
    oi_values = {i: set() for i in otherinfo}

    seen = set()
    rows = []
    func, exonic, clnsig = {}, {}, {}
    n_cadd = n_novel = n_spliceai = n_squirls = n_path = 0
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < len(key_idx):
            continue
        rows.append(parts)
        for i in otherinfo:
            if i < len(parts):
                oi_values[i].add(parts[i])
        key = tuple(parts[i] for i in key_idx)
        if key in seen:
            continue
        seen.add(key)

        def val(i):
            return parts[i] if i is not None and i < len(parts) else None

        if val(i_func) not in (None, "."):
            func[parts[i_func]] = func.get(parts[i_func], 0) + 1
        if val(i_exonic) not in (None, "."):
            exonic[parts[i_exonic]] = exonic.get(parts[i_exonic], 0) + 1
        v = val(i_clnsig)
        if v not in (None, "."):
            clnsig[v] = clnsig.get(v, 0) + 1
            if re.match(r"(Pathogenic|Likely_pathogenic|Pathogenic/Likely_pathogenic)($|\|)", v):
                n_path += 1
        if (c := numeric(val(i_cadd))) is not None and c >= CADD_MIN:
            n_cadd += 1
        if i_afs and all(val(i) in (".", None, "") for i in i_afs):
            n_novel += 1
        sai = [s for s in (numeric(val(i)) for i in i_spliceai) if s is not None]
        if sai and max(sai) >= SPLICE_MIN:
            n_spliceai += 1
        sq = [s for s in (numeric(val(i)) for i in i_squirls) if s is not None]
        if sq and max(sq) >= SPLICE_MIN:
            n_squirls += 1

    rare["rare_variants"] = len(seen)
    rare["rare_variant_sample_rows"] = len(rows)
    for label, counts in (("func", func), ("exonic_func", exonic)):
        for value in sorted(counts):
            rare[f"{label}:{value.replace(' ', '_')}"] = counts[value]

    for value in sorted(clnsig):
        flags[f"clnsig:{value.replace(' ', '_')}"] = clnsig[value]
    if i_clnsig is not None:
        flags["clinvar_path_or_likely_path"] = n_path
    if i_cadd is not None:
        flags[f"cadd_ge_{CADD_MIN:g}"] = n_cadd
    if i_afs:
        flags["gnomad_novel"] = n_novel
    if i_spliceai:
        flags[f"spliceai_ge_{SPLICE_MIN:g}"] = n_spliceai
    if i_squirls:
        flags[f"squirls_ge_{SPLICE_MIN:g}"] = n_squirls

    sample_cols = [i for i in otherinfo if oi_values[i] and oi_values[i] <= sample_names]
    if sample_cols:
        i_sample = sample_cols[-1]
        for parts in rows:
            if i_sample < len(parts):
                per_sample[parts[i_sample]] = per_sample.get(parts[i_sample], 0) + 1
    elif sample_names and otherinfo:
        note("no Otherinfo column holds only VCF sample names; "
             "per-sample rare counts are omitted")
    return rare, flags, per_sample


# --- somalier -----------------------------------------------------------------------------

def parse_samples_tsv(path):
    """{sample: (x_copies, y_copies, inferred_sex)} from X/Y depth ratios."""
    out = {}
    lines = read_lines(path, "somalier samples.tsv")
    if lines is None:
        return out
    i_s = header_index(lines[0], "sample_id")
    i_d = header_index(lines[0], "depth_mean")
    i_x = header_index(lines[0], "X_depth_mean")
    i_y = header_index(lines[0], "Y_depth_mean")
    if None in (i_s, i_d, i_x, i_y):
        note("samples.tsv lacks sample_id/depth_mean/X_depth_mean/Y_depth_mean; "
             "inferred sex is omitted")
        return out
    for line in lines[1:]:
        parts = line.split("\t")
        if not line.strip() or len(parts) <= max(i_s, i_d, i_x, i_y):
            continue
        depth = numeric(parts[i_d])
        if not depth:
            note(f"sample {parts[i_s]} has no usable depth_mean; its sex is undetermined")
            out[parts[i_s]] = (".", ".", "undetermined")
            continue
        xr = (numeric(parts[i_x]) or 0) / depth
        yr = (numeric(parts[i_y]) or 0) / depth
        x = 1 if xr <= X_HEMI_MAX else 2 if xr >= X_DIPLOID_MIN else "."
        y = 0 if yr <= Y_ABSENT_MAX else 1 if yr >= Y_PRESENT_MIN else "."
        sex = ("male" if (x, y) == (1, 1) else
               "female" if (x, y) == (2, 0) else "undetermined")
        out[parts[i_s]] = (x, y, sex)
    return out


def parse_ploidy(path):
    """{sample: (ped_sex, agreement)} from check_somalier's ploidy.tsv."""
    out = {}
    lines = read_lines(path, "ploidy.tsv")
    if lines is None:
        return out
    i_s = header_index(lines[0], "sample")
    i_p = header_index(lines[0], "ped_sex")
    i_a = header_index(lines[0], "agreement")
    if None in (i_s, i_p, i_a):
        note("ploidy.tsv lacks sample/ped_sex/agreement; the PED cross-check is omitted")
        return out
    for line in lines[1:]:
        parts = line.split("\t")
        if line.strip() and len(parts) > max(i_s, i_p, i_a):
            out[parts[i_s]] = (parts[i_p], parts[i_a])
    return out


def parse_ancestry(path):
    """{sample: (predicted, probability-or-None)} from the cohort ancestry TSV."""
    out = {}
    lines = read_lines(path, "somalier ancestry TSV")
    if lines is None:
        return out
    i_s = header_index(lines[0], "sample_id")
    i_p = header_index(lines[0], "predicted_ancestry")
    if None in (i_s, i_p):
        note("ancestry TSV lacks sample_id/predicted_ancestry; ancestry rows are omitted")
        return out
    for line in lines[1:]:
        parts = line.split("\t")
        if not line.strip() or len(parts) <= max(i_s, i_p):
            continue
        predicted = parts[i_p]
        i_prob = header_index(lines[0], f"{predicted}_prob")
        prob = parts[i_prob] if i_prob is not None and i_prob < len(parts) else None
        out[parts[i_s]] = (predicted, prob)
    return out


def parse_pairs(path):
    """(max_relatedness, "A,B") over pairs.tsv, or None."""
    lines = read_lines(path, "somalier pairs.tsv")
    if lines is None:
        return None
    i_a = header_index(lines[0], "sample_a")
    i_b = header_index(lines[0], "sample_b")
    i_r = header_index(lines[0], "relatedness")
    if None in (i_a, i_b, i_r):
        note("pairs.tsv lacks sample_a/sample_b/relatedness; max relatedness is omitted")
        return None
    best = None
    for line in lines[1:]:
        parts = line.split("\t")
        if not line.strip() or len(parts) <= max(i_a, i_b, i_r):
            continue
        rel = numeric(parts[i_r])
        if rel is not None and (best is None or rel > best[0]):
            best = (rel, f"{parts[i_a]},{parts[i_b]}")
    return best


# --- main ---------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stats", required=True, help="bcftools stats -s - output")
    ap.add_argument("--multianno", required=True, help="merged hg38_multianno.tsv")
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--data-type", required=True)
    ap.add_argument("--run-date", required=True)
    ap.add_argument("--clinvar-date", required=True)
    ap.add_argument("--samples-tsv", help="somalier relate samples.tsv")
    ap.add_argument("--ploidy", help="check_somalier ploidy.tsv")
    ap.add_argument("--ancestry", help="cohort somalier-ancestry.tsv")
    ap.add_argument("--pairs", help="somalier relate pairs.tsv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    sn, tstv, psc = parse_stats(args.stats)
    rare, flags, rare_per_sample = parse_multianno(args.multianno, set(psc))
    sexes = parse_samples_tsv(args.samples_tsv)
    ploidy = parse_ploidy(args.ploidy)
    ancestry = parse_ancestry(args.ancestry)
    max_rel = parse_pairs(args.pairs)

    rows = []

    def row(section, metric, sample, value):
        rows.append((section, metric, sample, str(value)))

    row("run", "cohort", ".", args.cohort)
    row("run", "data_type", ".", args.data_type)
    row("run", "run_date", ".", args.run_date)
    row("run", "clinvar_date", ".", args.clinvar_date)
    if "n_samples" in sn:
        row("run", "n_samples", ".", sn["n_samples"])

    for _key, metric in SN_METRICS:
        if metric != "n_samples" and metric in sn:
            row("cohort", metric, ".", sn[metric])
    if tstv is not None:
        row("cohort", "ts_tv", ".", tstv)

    for sample in sorted(psc):
        cols = psc[sample]
        for colname, metric in PSC_METRICS:
            if metric == "het_hom_ratio":
                het = numeric(cols.get("nHets"))
                hom = numeric(cols.get("nNonRefHom"))
                if het is not None and hom is not None:
                    row("sample", metric, sample,
                        f"{het / hom:.4f}" if hom else ".")
            elif colname in cols:
                row("sample", metric, sample, cols[colname])

    for metric in rare:  # insertion order: totals first, then sorted breakdowns
        row("rare_subset", metric, ".", rare[metric])
    for sample in sorted(rare_per_sample):
        row("rare_subset", "rare_variants", sample, rare_per_sample[sample])
    for metric in flags:
        row("flags", metric, ".", flags[metric])

    for sample in sorted(sexes):
        x, y, sex = sexes[sample]
        row("qc", "inferred_sex", sample, sex)
        row("qc", "x_copies", sample, x)
        row("qc", "y_copies", sample, y)
    for sample in sorted(ploidy):
        ped_sex, agreement = ploidy[sample]
        row("qc", "ped_sex", sample, ped_sex)
        row("qc", "sex_agreement", sample, agreement)
    for sample in sorted(ancestry):
        predicted, prob = ancestry[sample]
        row("qc", "predicted_ancestry", sample, predicted)
        if prob is not None:
            row("qc", "ancestry_prob", sample, prob)
    if max_rel is not None:
        row("qc", "max_relatedness", ".", f"{max_rel[0]:g}")
        row("qc", "max_relatedness_pair", ".", max_rel[1])

    with open(args.out, "w") as fh:
        fh.write("#section\tmetric\tsample\tvalue\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
