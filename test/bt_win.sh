#!/usr/bin/env bash
# bedtools-crossplatform: run a native Windows bedtools.exe from the MSYS2 bash
# test scripts (test.sh sets BT to this wrapper on Windows).
#
# A native program can't open bash process substitutions: <(...) becomes
# /proc/<pid>/fd/<n>, a path that only MSYS2 programs understand. This wrapper
# copies each such argument (inherited here as fd <n>) to a temporary file, maps
# /dev/stdin to bedtools' own "-", and passes everything else through unchanged.
BTEXE=${BTEXE:?set BTEXE to bedtools.exe}
tmps=()
args=()
for a in "$@"; do
    case $a in
        /proc/[0-9]*/fd/[0-9]*|/dev/fd/[0-9]*)
            t=$(mktemp)
            cat "/dev/fd/${a##*/}" > "$t"
            tmps+=("$t")
            args+=("$t")
            ;;
        /dev/stdin)
            args+=(-)
            ;;
        *)
            args+=("$a")
            ;;
    esac
done
"$BTEXE" "${args[@]}"
rc=$?
[ ${#tmps[@]} -eq 0 ] || rm -f "${tmps[@]}"
exit $rc
