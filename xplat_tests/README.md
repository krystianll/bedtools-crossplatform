# Cross-platform parity tests

`make test` checks each build against the upstream expected results. These
scripts add a byte-level check that a native Windows build writes exactly what
the macOS (or Linux) build writes for the same command lines.

- `parity.py dump BEDTOOLS_EXE DATA_DIR WORK_DIR OUT.json [case-substring]`
  runs every case in `CASES` with one build, each in its own empty directory
  under `WORK_DIR`, and records the exit code, stdout, stderr and every file the
  case wrote. `DATA_DIR` is the tree's `test/` directory (the upstream test data).
- `parity.py cmp A.json B.json` compares two dumps and prints `SAME`,
  `EXPECTED-DIFF` (listed in `EXPECTED_XPLAT`, with the reason) or `DIFF` with
  a diff per case, then the counts. The exit status is 1 if any case differs
  unexpectedly.
- `run_parity.sh BEDTOOLS_EXE OUT.json [WORK_DIR] [case-substring]` runs the
  dump with the tree's `test/` data and a temporary working directory.

What is compared:

- The working directory and `DATA_DIR` are replaced by `<CWD>` and `<DATA>`,
  in every spelling: as given and resolved, with `/` or `\`, and as `C:\…`,
  `C:/…` and `/c/…`. A `\` inside a path after these tags counts as `/`. All
  other bytes must be equal. Line endings are not normalised, so a CR the macOS
  build does not write is a difference (shown as `\r`).
- gzip and BGZF output (BAM) is compared after decompression, plus its type:
  the macOS build deflates with zlib, the Windows build with libdeflate, so the
  compressed bytes may differ. BAM is shown as header text, references and one
  hex line per record.
- The version (`v2.31.1-N-gHASH` from `git describe` in a checkout, or
  `v2.31.1` from an exported tree) is replaced by `v<VERSION>`.

The 186 cases cover every subcommand in `bedtools --help`, with the
common options of each (intersect `-wa/-wb/-v/-u/-c/-f/-r/-s/-split/-sorted`,
BAM in and out, CRAM; closest `-d/-D/-fu/-fd/-k`; coverage `-hist/-d/-mean`;
genomecov `-bg/-d/-ibam`; getfasta `-fo/-s/-name/-split/-tab`, a bgzipped FASTA,
a new `.fai`; seeded random/shuffle/sample), and these inputs: stdin (`-`,
`stdin`, gzip and BAM through stdin, empty stdin), gzip and bgzip files, CRLF
BED and FASTA, files without a final newline, empty files, missing files,
relative paths in a subdirectory, a directory with a space, a non-ASCII
directory (`zażółć`), absolute paths with `/` and with the platform's own
separator (file names in `multiinter -header`, `intersect -filenames`,
`links`), `--version`, `--help`, usage errors, and `bamtobed -i -` with empty
stdin.

## macOS / Linux

```bash
export CPPFLAGS="-I$(brew --prefix)/include" LDFLAGS="-L$(brew --prefix)/lib"  # macOS
make -j8
xplat_tests/run_parity.sh bin/bedtools mac.json
```

## Windows (MSYS2)

Build the tree as described in the main `README.md`, then in the CLANGARM64 or
UCRT64 shell (the build's DLL directory first in `PATH`, and that
environment's Python):

```bash
PATH=/clangarm64/bin:$PATH PYTHON=/clangarm64/bin/python \
    xplat_tests/run_parity.sh bin/bedtools.exe win.json
```

Copy `win.json` next to `mac.json` and compare:
`python3 xplat_tests/parity.py cmp win.json mac.json`.

## Expected differences (`EXPECTED_XPLAT`)

- `intersect-cram`: CRAM output. Its blocks are deflated by zlib or libdeflate
  and are not decoded here; `intersect-cram-bed` compares the same records as
  text.
- `split-size`: `bedtools split -a size` sorts the records by size with
  `std::sort`, which does not keep the order of equal sizes. libstdc++ (GCC: the
  UCRT64 build, and Linux) orders such ties differently from libc++ (clang:
  macOS and CLANGARM64), so some records land in other files. This is upstream
  behaviour, not a Windows problem.

Known open difference: the two non-ASCII path cases fail on Windows because
`bedtools.exe` gets its arguments in the ANSI code page (`zażółć` arrives as
`zazólc`). Embedding a UTF-8 code-page manifest
(`<activeCodePage>UTF-8</activeCodePage>`) in the executable would fix this.
