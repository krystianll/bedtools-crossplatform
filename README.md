# bedtools — patched, cross-platform build

A minimal patched copy of [**bedtools**](https://github.com/arq5x/bedtools2)
(upstream `master` at `614e9a5`, after release 2.31.1) that builds and runs
natively on **Linux, macOS (Intel + Apple Silicon), and Windows (x64 + ARM64)**,
and uses **htslib 1.24** instead of the htslib 1.9 bundled upstream.

For bedtools' documentation, options and canonical source, see the upstream
repository <https://github.com/arq5x/bedtools2> and the
[documentation](https://bedtools.readthedocs.io). **This README covers only how
to build this tree and what was changed**; everything else is unchanged from
upstream.

Upstream bedtools builds only on Unix-like systems. This tree adds a native
Windows build (MSYS2, no Cygwin runtime), fixes the Windows text-mode and 32-bit
`long` problems that broke its output there, and initialises several option
flags that upstream reads before setting (undefined behaviour on every
platform; visible on Windows as `closest`, `coverage` and `complement` acting
as if extra options were given).

The upstream test suite (`make test`) passes in full — 821 checks, only the
deliberate `negativecontrol` failing — on macOS, Linux aarch64, Windows ARM64
and Windows x64.

---

## Building from source

Run every command **from this directory** (the one containing `Makefile`). The
build makes one program, `bin/bedtools` (`bin/bedtools.exe` on Windows); every
tool is a subcommand (`bedtools intersect …`).

htslib (in `src/utils/htslib`) is compiled along with bedtools. Run its
`./configure` first to pick up optional libraries (libdeflate, …); without it,
htslib uses a default configuration that needs zlib, bzip2, liblzma and
libcurl.

### Linux

```bash
sudo apt install build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
make -j8
```

### macOS (Apple clang, Intel or Apple Silicon)

liblzma comes from Homebrew (`brew install xz`):

```bash
export CPPFLAGS="-I$(brew --prefix)/include" LDFLAGS="-L$(brew --prefix)/lib"
make -j8
```

### Windows x64 (MSYS2 **UCRT64** shell) and ARM64 (MSYS2 **CLANGARM64** shell)

Install the toolchain and libraries (`<env>` = `ucrt-x86_64` or
`clang-aarch64`):

```bash
pacman -S --needed make mingw-w64-<env>-gcc   # UCRT64; on CLANGARM64: mingw-w64-<env>-clang
pacman -S --needed mingw-w64-<env>-{zlib,bzip2,xz,libdeflate,libsystre,curl,openssl}
```

`libsystre` provides the POSIX `<regex.h>` htslib needs on Windows. Then:

```bash
# UCRT64
(cd src/utils/htslib && ./configure)
make -j8

# CLANGARM64
(cd src/utils/htslib && CC=clang ./configure)
make -j8 CC=clang CXX=clang++
```

The result is `bin/bedtools.exe`, a native Windows program; it needs the MSYS2
runtime DLLs it links (zlib1.dll, libcurl-4.dll, …) next to it or on `PATH`. The
other files in `bin/` (`intersectBed`, …) are upstream's one-line shell scripts
for the old command names; on Windows use `bedtools <subcommand>` instead.

### Verify (any platform)

```bash
./bin/bedtools --version
make test        # 821 "ok"; "Tools failing: negativecontrol" is expected
```

On Windows run `make test` from the same MSYS2 shell; `test.sh` runs
`bedtools.exe` through `test/bt_win.sh` (see the table below).

---

## Patches (vs upstream `614e9a5`)

Every edit to an upstream file carries a `bedtools-crossplatform:` comment, so

```bash
grep -rn "bedtools-crossplatform:" --exclude-dir=htslib src test Makefile
```

lists them. Windows-only changes are guarded by `#ifdef _WIN32` (or
`Windows_NT` in the Makefile and test scripts).

| File | Change | Reason |
|------|--------|--------|
| `src/utils/htslib/` | Replaced with the unmodified **htslib 1.24** release | Upstream bundles htslib 1.9 (2018). 1.24 has current CRAM/BGZF fixes and native ARM64 support, and matches the htslib of STAR-crossplatform and samtools 1.24, so an installer can ship one shared htslib. |
| `src/utils/BamTools/HtsCallback.{h,cpp}` *(new)*, `BamTools/Makefile.frag`, `SamHeader.hpp` | The stream-callback I/O (`hopen_callback`, `hts_open_callback`) that bedtools had added to its htslib 1.9, re-implemented on htslib's exported backend API (`hfile_init`, `hts_hopen`); a NULL check added | Keeps htslib unmodified. Used by `BamReader.hpp`, `SamHeader.hpp` and `InputStreamMgr.cpp`. |
| `Makefile` | Link the libraries htslib reports in `htslib_static.mk`; on Windows also `-lregex` | Follows htslib's configure result (libdeflate, regex, Winsock, OpenSSL on Windows). htslib's configure needs `-lregex` on MinGW but leaves it out of `htslib_static.mk`. |
| `src/bedtools.cpp` | `#ifdef _WIN32`: at the top of `main()`, `_fmode = _O_BINARY` and binary stdin/stdout/stderr | Windows text mode turned every `\n` into `\r\n` and could alter BAM/BGZF bytes before detection. `_fmode` makes every file opened afterwards binary too (`getfasta -fo`, `split`, `bamtofastq -fq`, all inputs), so no call site needed changing. |
| `src/utils/Contexts/Context{Base,Closest,Complement,Coverage,Intersect}.cpp`, `src/utils/bedFile/bedFile.cpp` | Initialise `_runToQueryEnd`, `_isCram`, `_forceUpstream`, `_forceDownstream`, `_onlyChromsWithBedRecords`, `_mean`, `_shouldRunToDbEnd`, `_lineNum` | Read before being set (undefined behaviour everywhere). On Windows plain `closest` acted as `-fu`/`-fd`, `coverage` as `-mean`, `complement` as `-L`, and `intersect` could drain DB records before the sorted-input check. |
| `src/utils/GenomeFile/GenomeFile.cpp` | Start `_genomeLength` at 0; the BAM-header constructor builds offsets and length like the file constructor; parse chromosome sizes with `strtoll` | The uninitialised total broke `shuffle`/`sample`; `long` is 32 bits on Windows, so `atol` truncated large genomes. |
| `src/coverageFile/coverageFile.cpp` | `#ifdef _WIN32`: a local `asprintf` | MinGW has none. |
| `src/fisher/kfunc.cpp` | `M_SQRT2` fallback | MinGW's `math.h` hides it in strict C++ mode. |
| `src/utils/Fasta/Fasta.h` | `<sys/mman.h>` only off Windows | Absent on Windows; `Fasta.cpp` no longer memory-maps. |
| `src/utils/Fasta/LargeFileSupport.h`, `src/nucBed/LargeFileSupport.h` | `int64_t` instead of `__int64_t` in the `WIN32` branch | `__int64_t` is not a MinGW type. |
| `src/regressTest/RegressTest.cpp` | `_mkdir` on Windows | Windows' `mkdir` takes no mode. |
| `sequenceUtils.h`, `lineFileUtilities.h`, `GenomeFile.h`, `SingleLineDelimTextFileReader.h`, `BamTools/include/api/BamAux.h` | `#include <stdint.h>` | They use `int64_t`/`int32_t` without it; gcc 16 does not get it indirectly. |
| `test/test.sh` | On Windows, set `BT` to `test/bt_win.sh`; use `python3` where there is no `python` | Lets the unchanged test scripts drive `bedtools.exe`; some tests run `python`. |
| `test/bt_win.sh` *(new)* | Wrapper that copies bash process substitutions (`<(...)`, which appear as `/proc/<pid>/fd/<n>`) to temporary files and maps `/dev/stdin` to `-` | A native Windows program cannot open those MSYS2-only paths. |
| `test/fisher/test-fisher.sh` | `fisher.t5`: drop the backslashes from the long file name on Windows | `\` is a path separator there. |
| `xplat_tests/` *(new)* | Parity harness: dumps 186 cases (every subcommand; stdin, gzip/bgzip, CRLF input, spaces, non-ASCII and backslash paths) and compares a Windows dump with a macOS one byte for byte, BAM/BGZF after decompression | Checks that the Windows build behaves exactly like macOS, beyond the pass/fail of `make test`. See `xplat_tests/README.md`. |

### Upstream behaviour kept as it is

- **CRLF input is not normalised.** Files are read in binary mode, as on Unix, so
  a BED/GFF/VCF file with Windows line endings keeps a trailing `\r` in its last
  field — exactly what upstream bedtools does on Linux. Convert such files to LF
  first (`dos2unix`).
- **Non-ASCII paths on Windows** (outside the ANSI code page) need a UTF-8
  application manifest (`<activeCodePage>UTF-8</activeCodePage>`) embedded in
  `bedtools.exe`, as it uses the narrow-character file APIs; without it
  `zażółć/ą.bed` reaches bedtools as `zazólc/a.bed`. The ngs-tools installer
  embeds one.
- **`split -a size`** orders equal-sized records with an unstable sort, so the
  files they land in can differ between C++ libraries (libstdc++ vs libc++).
- **`bedtools regresstest`** (a developer tool) uses POSIX shell commands and does
  not work on Windows.
- **HTTPS/S3 input** goes through htslib's libcurl. On Windows, point
  `CURL_CA_BUNDLE` at a CA certificate file (the ngs-tools installer does) or
  certificate checks fail.

---

## Related work

[win-ngs](https://github.com/win-ngs) publishes community Windows builds of
several NGS tools, including
[bedtools-windows-build](https://github.com/win-ngs/bedtools-windows-build)
(bedtools 2.31.1, x64, MSYS2 UCRT64). Its README's list of Windows problems —
text-mode I/O, the 32-bit genome sizes, and the uninitialised option flags it
traced through the upstream test suite — was a useful checklist for this fork;
this fork shares no code with it. It differs in building natively for ARM64 as
well as x64, using htslib 1.24, using one `_fmode` switch instead of per-call
binary modes, and running the upstream test scripts unchanged through a small
wrapper.

## Redistribution & licensing

bedtools is **MIT** licensed (see `LICENSE`); htslib is MIT/BSD
(`src/utils/htslib/LICENSE`). A distributed binary also carries the notices of
the libraries linked in — zlib, bzip2, liblzma, libdeflate, libcurl and its TLS
stack, and on Windows the libsystre/TRE/gettext/libiconv regex chain and the
GCC/libstdc++ (or LLVM libc++) and winpthread runtimes. Ship each DLL's license
with it.
