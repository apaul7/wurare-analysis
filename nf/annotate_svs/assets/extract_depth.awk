# Extract duphold's per-sample depth fields into a flat table keyed on record ID.
#
# Usage:  awk -f extract_depth.awk single_sample_duphold.vcf
# Emits:  ID <tab> SAMPLE <tab> DHFFC <tab> DHBFC <tab> DHFC <tab> DHBZ
#
# Written as awk rather than `bcftools query -f '[%DHFFC...]'` because bcftools hard-errors
# on a tag the header does not declare -- "no such tag defined in the VCF header:
# FORMAT/DHBZ" -- and duphold does not always emit all four. DHBZ in particular needs enough
# GC-matched bins to exist, so it is absent on small references and on some real ones.
# Reading the FORMAT column by name is tolerant by construction: a tag that is not there
# becomes ".", which is what the recombination wants anyway.
#
# The sample name is taken from the VCF's own #CHROM line, never from the filename.

BEGIN { FS = OFS = "\t" }

/^#CHROM/ { sample = $10; next }
/^#/ { next }

{
    n = split($9, keys, ":")
    nv = split($10, vals, ":")
    delete f
    # VCF 4.2 section 1.6.2 lets a caller drop trailing FORMAT subfields, so the value count can be
    # shorter than the key count. Reading past the end yields the empty string, not ".", and
    # an empty subfield is invalid VCF once merge_depth writes it back into a sample column.
    # Short means absent, which is what an undeclared tag already means here.
    for (i = 1; i <= n; i++) f[keys[i]] = (i <= nv && vals[i] != "") ? vals[i] : "."

    print $3, sample,
          (("DHFFC" in f) ? f["DHFFC"] : "."),
          (("DHBFC" in f) ? f["DHBFC"] : "."),
          (("DHFC"  in f) ? f["DHFC"]  : "."),
          (("DHBZ"  in f) ? f["DHBZ"]  : ".")
}
