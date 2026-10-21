#!/usr/bin/env bash
# Dump the results of every parity case for one bedtools build (see README.md).
# usage: xplat_tests/run_parity.sh BEDTOOLS_EXE OUT.json [WORK_DIR] [case-substring]
#   BEDTOOLS_EXE  the build to run (bin/bedtools, or bin/bedtools.exe on Windows)
#   OUT.json      where to write the dump; compare two dumps with
#                 python3 xplat_tests/parity.py cmp A.json B.json
#   WORK_DIR      scratch directory (default: a new temporary one, removed after)
# The data are the tree's test/ directory. PYTHON overrides the interpreter
# (default python3, else python). On Windows (MSYS2) put the build's DLL
# directory first in PATH, e.g. PATH=/clangarm64/bin:$PATH.
set -eu
[ $# -ge 2 ] || { sed -n '2,10p' "$0"; exit 2; }
H=$(cd "$(dirname "$0")/.." && pwd)
EXE=$1 OUT=$2 W=${3:-} ONLY=${4:-}
PY=${PYTHON:-$(command -v python3 || command -v python)}
[ -x "$EXE" ] || { echo "not executable: $EXE"; exit 2; }
if [ -z "$W" ]; then
	W=$(mktemp -d)
	trap 'rm -rf "$W"' EXIT
fi
"$PY" "$H/xplat_tests/parity.py" dump "$EXE" "$H/test" "$W" "$OUT" ${ONLY:+"$ONLY"}
