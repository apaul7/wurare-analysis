nextflow.enable.types = true

// Parsing for the two input sheets.
//
// Sample IDs are deliberately absent from the VCF sheet: they are read from each VCF
// header at preflight (Phase 1), because a sheet that restates them is a second source of
// truth that will drift, and drift here means one sample silently becoming two columns.

record SvInput {
    sample_set: String
    caller: String
    joint: Boolean
    vcf: Path
    tbi: Path
}

record AlignmentEntry {
    sample: String
    alignment: Path
    index: Path
}

// The caller becomes the SVDB `--vcf file:tag` tag, and the tag and priority lists are
// colon- and comma-delimited, so a caller containing either corrupts the merge invocation
// rather than failing.
//
// Allow-list rather than deny-list: both values reach a shell. `caller` is interpolated
// UNQUOTED into svdb's --vcf/--priority arguments and `sample_set` into output filenames, so
// "manta v2" alone splits into two svdb arguments. A trailing space from a spreadsheet export
// is the likeliest way to meet this.
def check_sheet_token(String value, String field, String sample_set) {
    if (!(value ==~ /[A-Za-z0-9._-]+/)) {
        error "--vcfs sheet: row '${sample_set}' has ${field}='${value}'; only letters, " +
              "digits, '.', '_' and '-' are allowed -- this value reaches both a shell " +
              "command line and an output filename"
    }
}

// Not a style preference: two rows sharing a (sample_set, caller) pair stage two files with
// the SAME name, and Nextflow keeps one symlink and silently drops the rest (measured in this
// repo -- see nf/shared/modules/versions.nf). One callset is then annotated, no warning.
def check_unique_rows(List rows) {
    def seen = [:]
    rows.each { r ->
        def key = "${r.sample_set}\u0000${r.caller}"
        if (seen.containsKey(key)) {
            error "--vcfs sheet: duplicate row for sample_set='${r.sample_set}' " +
                  "caller='${r.caller}'. Both would be staged under the same filename and " +
                  "one would be silently dropped; give the second a distinct caller name"
        }
        seen[key] = true
    }
    return rows
}

// Checked by name rather than position so a reordered or renamed sheet fails naming the
// column, instead of silently reading `caller` out of the `joint` slot.
def require_columns(Map row, List<String> columns, String sheet) {
    def missing = columns.findAll { c -> !row.containsKey(c) || !row[c] }
    if (missing) {
        error "${sheet}: row is missing required value(s) ${missing.join(', ')} -- got ${row}"
    }
}

// `joint` cannot be inferred from sample count: a multi-sample VCF produced by
// concatenating single-sample calls is not joint-called, and treating it as such would
// give it priority over the callset whose genotypes are actually worth protecting.
def parse_joint(String value, String sample_set) {
    def v = value.trim().toLowerCase()
    if (v in ['true', 'yes', '1'])  return true
    if (v in ['false', 'no', '0'])  return false
    error "--vcfs sheet: row '${sample_set}' has joint='${value}'; expected true or false"
}

def to_sv_input(Map row) {
    require_columns(row, ['sample_set', 'caller', 'joint', 'vcf', 'tbi'], '--vcfs sheet')
    check_sheet_token((row.caller as String).trim(), 'caller',
                      (row.sample_set as String).trim())
    check_sheet_token((row.sample_set as String).trim(), 'sample_set',
                      (row.sample_set as String).trim())
    new SvInput(
        sample_set: (row.sample_set as String).trim(),
        caller:     (row.caller as String).trim(),
        joint:      parse_joint(row.joint as String, row.sample_set as String),
        vcf:        file((row.vcf as String).trim(), checkIfExists: true),
        tbi:        file((row.tbi as String).trim(), checkIfExists: true)
    )
}

def to_alignment(Map row) {
    require_columns(row, ['sample', 'alignment', 'alignment_index'], '--alignments sheet')
    new AlignmentEntry(
        sample:    (row.sample as String).trim(),
        alignment: file((row.alignment as String).trim(), checkIfExists: true),
        index:     file((row.alignment_index as String).trim(), checkIfExists: true)
    )
}
