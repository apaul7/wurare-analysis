#!/usr/bin/env bash
#
# Runs every tests/test_*.py and summarises. Each test is a standalone script that decides for
# itself whether it can run: the ones needing `nextflow` or `pysam` print "SKIP: ..." and exit
# 0 rather than failing, so a laptop without the full toolchain still gets a useful signal.
#
#   ./tests/run.sh                  # everything
#   ./tests/run.sh test_add_scores  # one, by name (with or without the test_ prefix or .py)
#
# Exit status is the number of FAILED tests, so `./tests/run.sh && echo green` works.
#
# Existed as nothing at all before -- tests were run one file at a time by hand, which meant
# "do they all still pass" had no answer short of twenty invocations.

set -uo pipefail

cd "$(dirname "$0")/.."
PY=${PYTHON:-python3}

if [ $# -gt 0 ]; then
    targets=()
    for a in "$@"; do
        a=${a%.py}; a=${a#test_}
        targets+=("tests/test_${a}.py")
    done
else
    targets=(tests/test_*.py)
fi

# EXCLUDE="test_foo test_bar" drops tests by name -- CI uses this to keep the
# docker e2e tier out of the fast job, where it would otherwise pull gigabytes.
if [ -n "${EXCLUDE:-}" ]; then
    kept=()
    for t in "${targets[@]}"; do
        case " $EXCLUDE " in
            *" $(basename "$t" .py) "*) ;;
            *) kept+=("$t") ;;
        esac
    done
    targets=("${kept[@]}")
fi

pass=0; fail=0; skip=0
failed_names=()
log=$(mktemp)
trap 'rm -f "$log"' EXIT

for t in "${targets[@]}"; do
    name=$(basename "$t" .py)
    if [ ! -f "$t" ]; then
        printf 'FAIL %s  (no such test)\n' "$name"
        fail=$((fail+1)); failed_names+=("$name")
        continue
    fi
    if "$PY" "$t" > "$log" 2>&1; then
        # A test that could not run says so on its first line rather than reporting a pass it
        # did not earn. Counted separately: a suite that is entirely skips is not a green one.
        if head -n1 "$log" | grep -q '^SKIP:'; then
            reason=$(head -n1 "$log" | sed 's/^SKIP: //')
            printf 'skip %s  -- %s\n' "$name" "$reason"
            # Under GitHub Actions a skip is yellow, not silently green: each one
            # becomes a warning annotation on the run.
            [ -n "${GITHUB_ACTIONS:-}" ] && printf '::warning title=skipped %s::%s\n' "$name" "$reason"
            skip=$((skip+1))
        else
            printf 'ok   %s\n' "$name"
            pass=$((pass+1))
        fi
    else
        printf 'FAIL %s\n' "$name"
        # The tests print their "FAIL <check>" verdicts before dumping tool output,
        # so a bare tail can scroll the actual reason out of view. Show both.
        grep -a '^ *FAIL' "$log" | head -10 | sed 's/^/       | /'
        sed 's/^/       | /' "$log" | tail -30
        fail=$((fail+1)); failed_names+=("$name")
    fi
done

printf '\n%d passed, %d failed, %d skipped\n' "$pass" "$fail" "$skip"
[ "$fail" -gt 0 ] && printf 'failed: %s\n' "${failed_names[*]}"
exit "$fail"
