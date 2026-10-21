set -e; # Alert user to any uncaught error

ulimit -c unlimited

STARTWD=$(pwd);

# bedtools-crossplatform: on Windows, run the native bedtools.exe through
# bt_win.sh, which turns bash process substitutions and /dev/stdin into paths a
# native program can open.
if [ "$OS" = "Windows_NT" ] && [ -z "${BT+x}" ]; then
    export BTEXE="${STARTWD}/../bin/bedtools.exe"
    export BT="${STARTWD}/bt_win.sh"
fi

# bedtools-crossplatform: some tests run "python"; use python3 where only that exists.
if ! command -v python >/dev/null 2>&1 && command -v python3 >/dev/null 2>&1; then
    PYSHIM=$(mktemp -d)
    ln -s "$(command -v python3)" "$PYSHIM/python"
    export PATH="$PYSHIM:$PATH"
fi
TOOL_PASSES="";
TOOL_FAILURES="";

for tool in $(ls); do
    [ -d "${STARTWD}/${tool}" ] || continue;
    echo "Testing bedtools $tool:";
    cd "${STARTWD}/${tool}";
    bash "test-${tool}.sh" \
        && TOOL_PASSES="$TOOL_PASSES $tool" \
        || TOOL_FAILURES="$TOOL_FAILURES $tool";
done

echo
echo
echo "--------------------------"
echo " Test Results             "
echo "--------------------------"
echo "Tools passing: $TOOL_PASSES"
echo "Tools failing: $TOOL_FAILURES"
echo "NB: the 'negativecontrol' test is supposed to fail. If it wasn't caught, "
echo "something went wrong with this test script."
[ "$TOOL_FAILURES" = " negativecontrol" ] || exit 1;
