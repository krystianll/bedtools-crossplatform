"""Cross-platform parity of bedtools: run the same command lines with two builds
(e.g. native Windows and macOS) and compare the results byte for byte.

usage: python3 parity.py dump BEDTOOLS_EXE DATA_DIR WORK_DIR OUT.json [case-substring]
       python3 parity.py cmp A.json B.json      (e.g. Windows vs macOS)

DATA_DIR is the tree's test/ directory (the upstream test data). `dump` runs
every case in its own empty directory under WORK_DIR and records exit code,
stdout, stderr and every file the case wrote there. The working directory and
DATA_DIR are replaced by <CWD> and <DATA> (as given or resolved, with either
separator, and in the C:\\ C:/ /c/ drive forms). Line endings are NOT
normalised: a CR that macOS does not write is a difference. gzip/BGZF output is
recorded as its compression type plus the decompressed bytes (the builds use
different deflate libraries); BAM additionally as header text, references and
one hex line per record. `cmp` prints SAME, EXPECTED-DIFF (listed in
EXPECTED_XPLAT with the reason) or DIFF per case; exit status 1 means an
unexpected difference.
"""
import difflib, gzip, json, os, re, shutil, subprocess, sys

D, C = "{D}", "{CWD}"
NZ = C + "/za\u017c\u00f3\u0142\u0107"          # a non-ASCII directory
SP = C + "/a b"                                 # a directory with a space
I = D + "/intersect/"
G19 = D + "/sample/human.hg19.genome"
TINY = D + "/slop/tiny.genome"
FA = D + "/getfasta/t.fa"

# Inputs generated in the working directory (bytes, LF unless stated).
BEDPE = (b"chr1\t5\t25\tchr1\t90\t120\tp1\t30\t+\t-\n"
         b"chr1\t300\t400\tchr1\t500\t600\tp2\t10\t+\t+\n"
         b"chr1\t1000\t1100\tchr2\t10\t20\tp3\t5\t-\t+\n")
BEDPE_B = (b"chr1\t0\t30\tchr1\t100\t130\tq1\t1\t+\t-\n"
           b"chr1\t250\t350\tchr1\t550\t650\tq2\t2\t+\t+\n")
BG = [b"chr1\t0\t100\t1\nchr1\t100\t200\t2\nchr2\t0\t50\t3\n",
      b"chr1\t50\t150\t5\nchr2\t20\t80\t1.5\n",
      b"chr1\t0\t300\t0.25\n"]
OVL = b"chr1\t10\t20\tchr1\t15\t30\nchr1\t10\t20\tchr1\t30\t40\nchr2\t0\t100\tchr2\t50\t60\n"
GENOME2 = b"chr1\t2000\nchr2\t500\n"
REGIONS = b"chr1\t1\t10\tr1\t0\t+\nchr1\t5\t30\tr2\t0\t-\nchr1\t40\t50\tr3\t0\t+\n"   # for getfasta/t.fa
R = {"r.bed": REGIONS}


def crlf(path):
    """The file with every LF turned into CRLF (setup spec)."""
    return ("crlf", path)


def gz(path):
    return ("gzip", path)


# (name, args, options). Options: stdin = DATA-relative path, bytes, or b"" for
# empty input; setup = {working-dir-relative destination: source}, the source a
# {D} path (copied), bytes, crlf(path) or gz(path); env = extra environment;
# native = on Windows pass the {D}/{CWD} paths with backslashes only (otherwise
# they are written with /, mixed with the separators of DATA_DIR/WORK_DIR).
# Every argument and value may use {D} (DATA_DIR) and {CWD}.
CASES = [
    # the front end
    ("noargs", [], {}),
    ("help", ["--help"], {}),
    ("version", ["--version"], {}),
    ("contact", ["--contact"], {}),
    ("unknown-subcommand", ["nosuchtool"], {}),
    ("intersect-noargs", ["intersect"], {}),
    ("intersect-help", ["intersect", "-h"], {}),
    # intersect
    ("intersect", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed"], {}),
    ("intersect-wa-wb", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-wa", "-wb"], {}),
    ("intersect-wo", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-wo"], {}),
    ("intersect-v", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-v"], {}),
    ("intersect-u", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-u"], {}),
    ("intersect-c", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-c"], {}),
    ("intersect-f-r", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-f", "0.5", "-r", "-wao"], {}),
    ("intersect-s", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-s"], {}),
    ("intersect-S", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-S"], {}),
    ("intersect-split", ["intersect", "-a", I + "blocks.bed12", "-b", I + "a.bed", "-split", "-wa", "-wb"], {}),
    ("intersect-sorted", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-sorted"], {}),
    ("intersect-sorted-unsorted", ["intersect", "-a", I + "recordsOutOfOrder.bed", "-b", I + "a.bed", "-sorted"], {}),
    ("intersect-multi-b-names", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", I + "c.bed",
                                 "-names", "B", "C", "-wa", "-wb"], {}),
    ("intersect-multi-b-filenames", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", I + "c.bed",
                                     "-filenames", "-wa", "-wb"], {}),
    ("intersect-header", ["intersect", "-a", D + "/general/a.trackheader.bed",
                          "-b", I + "b.bed", "-header"], {}),
    ("intersect-abam", ["intersect", "-abam", I + "a.bam", "-b", I + "b.bed"], {}),
    ("intersect-abam-ubam", ["intersect", "-abam", I + "a.bam", "-b", I + "b.bed", "-ubam"], {}),
    ("intersect-abam-bed", ["intersect", "-abam", I + "a.bam", "-b", I + "b.bed", "-bed", "-wa", "-wb"], {}),
    ("intersect-bam-a", ["intersect", "-a", I + "a.bam", "-b", I + "b.bed", "-v"], {}),
    ("intersect-cram", ["intersect", "-a", I + "a.cram", "-b", I + "b.cram"],
     {"env": {"CRAM_REFERENCE": I + "test_ref.fa"}}),
    ("intersect-cram-bed", ["intersect", "-a", I + "a.cram", "-b", I + "b.cram", "-bed"],
     {"env": {"CRAM_REFERENCE": I + "test_ref.fa"}}),
    ("intersect-vcf", ["intersect", "-a", I + "a_vcfSVtest.vcf", "-b", I + "b_vcfSVtest.vcf"], {}),
    ("intersect-gff-vcf", ["intersect", "-a", I + "a.issue311.bed", "-b", I + "b.issue311.vcf", "-wa", "-wb"], {}),
    ("intersect-gzip", ["intersect", "-a", I + "a_gzipped.bed.gz", "-b", I + "b.bed"], {}),
    ("intersect-bgzip", ["intersect", "-a", I + "a_bgzipped.bed.gz", "-b", I + "b.bed"], {}),
    ("intersect-gzip-header", ["intersect", "-a", I + "a_withLargeHeader_gzipped.bed.gz", "-b", I + "b.bed",
                               "-header"], {}),
    ("intersect-gzip-made", ["intersect", "-a", "a.bed.gz", "-b", I + "b.bed", "-wa"],
     {"setup": {"a.bed.gz": gz(I + "a.bed")}}),
    ("intersect-stdin", ["intersect", "-a", "-", "-b", I + "b.bed"], {"stdin": "intersect/a.bed"}),
    ("intersect-stdin-word", ["intersect", "-a", "stdin", "-b", I + "b.bed"], {"stdin": "intersect/a.bed"}),
    ("intersect-stdin-gzip", ["intersect", "-a", "-", "-b", I + "b.bed"], {"stdin": "intersect/a_gzipped.bed.gz"}),
    ("intersect-stdin-bam", ["intersect", "-abam", "-", "-b", I + "b.bed", "-bed"], {"stdin": "intersect/a.bam"}),
    ("intersect-crlf-made", ["intersect", "-a", "crlf.bed", "-b", I + "b.bed", "-wa", "-wb"],
     {"setup": {"crlf.bed": crlf(I + "a.bed")}}),
    ("intersect-crlf-b", ["intersect", "-a", I + "a.bed", "-b", "crlf.bed", "-wa", "-wb"],
     {"setup": {"crlf.bed": crlf(I + "b.bed")}}),
    ("intersect-crlf-upstream", ["intersect", "-a", I + "dosLineChar_a.bed", "-b", I + "dosLineCharWithExtraTab_b.bed"], {}),
    ("intersect-no-final-newline", ["intersect", "-a", I + "oneRecordNoNewline.bed", "-b", I + "a.bed"], {}),
    ("intersect-empty-a", ["intersect", "-a", D + "/general/empty.bed", "-b", I + "b.bed"], {}),
    ("intersect-empty-stdin", ["intersect", "-a", "-", "-b", I + "b.bed"], {"stdin": b""}),
    ("intersect-missing-a", ["intersect", "-a", "nosuch.bed", "-b", I + "b.bed"], {}),
    ("intersect-missing-abs", ["intersect", "-a", C + "/nodir/nosuch.bed", "-b", I + "b.bed"], {}),
    ("intersect-relative-subdir", ["intersect", "-a", "sub/a.bed", "-b", "sub/b.bed", "-filenames", "-wa", "-wb"],
     {"setup": {"sub/a.bed": I + "a.bed", "sub/b.bed": I + "b.bed"}}),
    ("intersect-space-path", ["intersect", "-a", SP + "/x y.bed", "-b", SP + "/b b.bed", "-filenames", "-wa", "-wb"],
     {"setup": {"a b/x y.bed": I + "a.bed", "a b/b b.bed": I + "b.bed"}}),
    ("intersect-nonascii-path", ["intersect", "-a", NZ + "/\u0105.bed", "-b", NZ + "/\u0119.bed", "-filenames",
                                 "-wa", "-wb"],
     {"setup": {"za\u017c\u00f3\u0142\u0107/\u0105.bed": I + "a.bed", "za\u017c\u00f3\u0142\u0107/\u0119.bed": I + "b.bed"}}),
    ("intersect-filenames-native", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", SP + "/c.bed", "-filenames",
                                    "-wa", "-wb"], {"setup": {"a b/c.bed": I + "c.bed"}, "native": True}),
    ("intersect-missing-native", ["intersect", "-a", C + "/nodir/nosuch.bed", "-b", I + "b.bed"], {"native": True}),
    ("intersect-g-sorted", ["intersect", "-a", I + "a.bed", "-b", I + "b.bed", "-sorted", "-g", I + "human.hg19.genome"], {}),
    # window, closest, subtract
    ("window", ["window", "-a", I + "a.bed", "-b", I + "b.bed", "-w", "100"], {}),
    ("window-l-r-sm", ["window", "-a", I + "a.bed", "-b", I + "b.bed", "-l", "50", "-r", "10", "-sm"], {}),
    ("window-abam", ["window", "-abam", I + "a.bam", "-b", I + "b.bed", "-bed"], {}),
    ("closest", ["closest", "-a", D + "/closest/a.bed", "-b", D + "/closest/b.bed"], {}),
    ("closest-d", ["closest", "-a", D + "/closest/a.bed", "-b", D + "/closest/b.bed", "-d"], {}),
    ("closest-D-fu", ["closest", "-a", D + "/closest/strand-test-a.bed", "-b", D + "/closest/strand-test-b.bed",
                      "-D", "a", "-fu"], {}),
    ("closest-D-fd", ["closest", "-a", D + "/closest/strand-test-a.bed", "-b", D + "/closest/strand-test-b.bed",
                      "-D", "ref", "-fd"], {}),
    ("closest-k", ["closest", "-a", D + "/closest/close-a.bed", "-b", D + "/closest/close-b.bed", "-k", "3", "-d"], {}),
    ("subtract", ["subtract", "-a", D + "/subtract/a.bed", "-b", D + "/subtract/b.bed"], {}),
    ("subtract-A", ["subtract", "-a", D + "/subtract/a.bed", "-b", D + "/subtract/b.bed", "-A"], {}),
    # coverage, map, genomecov
    ("coverage", ["coverage", "-a", D + "/coverage/a.bed", "-b", D + "/coverage/b.bed"], {}),
    ("coverage-hist", ["coverage", "-a", D + "/coverage/a.bed", "-b", D + "/coverage/b.bed", "-hist"], {}),
    ("coverage-d", ["coverage", "-a", D + "/coverage/a.bed", "-b", D + "/coverage/b.bed", "-d"], {}),
    ("coverage-mean", ["coverage", "-a", D + "/coverage/a.bed", "-b", D + "/coverage/b.bed", "-mean"], {}),
    ("coverage-bam-split", ["coverage", "-a", D + "/coverage/three_blocks_match.bed", "-b",
                            D + "/coverage/three_blocks_match.bam", "-split"], {}),
    ("map", ["map", "-a", D + "/map/ivls.bed", "-b", D + "/map/values.bed", "-c", "5", "-o", "mean"], {}),
    ("map-multi-ops", ["map", "-a", D + "/map/ivls.bed", "-b", D + "/map/values.bed", "-c", "5,5,4,5",
                       "-o", "sum,stdev,collapse,median"], {}),
    ("map-bam", ["map", "-a", D + "/map/d.bed", "-b", D + "/map/fullFields.bam", "-c", "5", "-o", "mean"], {}),
    ("map-vcf-header", ["map", "-a", D + "/map/ivls.bed", "-b", D + "/map/test.vcf", "-header", "-c", "10",
                        "-o", "collapse"], {}),
    ("genomecov", ["genomecov", "-i", D + "/genomecov/y.bed", "-g", D + "/genomecov/genome.txt"], {}),
    ("genomecov-bg", ["genomecov", "-i", D + "/genomecov/y.bed", "-g", D + "/genomecov/genome.txt", "-bg"], {}),
    ("genomecov-d", ["genomecov", "-i", D + "/genomecov/y.bed", "-g", D + "/genomecov/genome.txt", "-d"], {}),
    ("genomecov-ibam", ["genomecov", "-ibam", D + "/genomecov/y.bam"], {}),
    ("genomecov-ibam-bg-split", ["genomecov", "-ibam", D + "/genomecov/y.bam", "-bg", "-split", "-scale", "0.5"], {}),
    ("genomecov-ibam-empty-cram", ["genomecov", "-ibam", D + "/genomecov/empty.cram"],
     {"env": {"CRAM_REFERENCE": D + "/genomecov/test_ref.fa"}}),
    # merge, cluster, complement, sort, spacing
    ("merge", ["merge", "-i", D + "/merge/a.bed"], {}),
    ("merge-d-c-o", ["merge", "-i", D + "/merge/a.names.bed", "-d", "10", "-c", "4,4", "-o", "collapse,count"], {}),
    ("merge-s", ["merge", "-i", D + "/merge/mixedStrands.bed", "-s", "-c", "6", "-o", "distinct"], {}),
    ("merge-precision", ["merge", "-i", D + "/merge/precisionTest.bed", "-c", "5", "-o", "mean"], {}),
    ("merge-bam", ["merge", "-i", D + "/merge/fullFields.bam", "-c", "5", "-o", "mean"], {}),
    ("merge-unsorted", ["merge", "-i", D + "/merge/unsorted.bed"], {}),
    ("merge-gzip-chained", ["merge", "-i", D + "/merge/chained.bed.gz"], {}),
    ("merge-empty", ["merge", "-i", D + "/general/empty.bed"], {}),
    ("merge-crlf", ["merge", "-i", "crlf.bed", "-c", "4", "-o", "collapse"],
     {"setup": {"crlf.bed": crlf(D + "/merge/a.names.bed")}}),
    ("cluster", ["cluster", "-i", D + "/cluster/in.bed"], {}),
    ("cluster-d-s", ["cluster", "-i", D + "/cluster/in.bed", "-d", "1000", "-s"], {}),
    ("complement", ["complement", "-i", D + "/complement/issue_503.bed", "-g", D + "/complement/issue_503.genome"], {}),
    ("sort", ["sort", "-i", D + "/sort/a.bed"], {}),
    ("sort-sizeD", ["sort", "-i", D + "/sort/a.bed", "-sizeD"], {}),
    ("sort-chrThenScoreA", ["sort", "-i", D + "/sort/b.bed", "-chrThenScoreA"], {}),
    ("sort-faidx", ["sort", "-i", D + "/sort/a.bed", "-faidx", D + "/sort/names.txt"], {}),
    ("sort-header-stdin", ["sort", "-i", "-", "-header"], {"stdin": "general/a.trackheader.bed"}),
    ("spacing", ["spacing", "-i", D + "/spacing/a.bed"], {}),
    # interval arithmetic with a genome file
    ("slop", ["slop", "-i", D + "/slop/a.bed", "-g", TINY, "-b", "5"], {}),
    ("slop-pct-s", ["slop", "-i", D + "/slop/a.bed", "-g", TINY, "-l", "0.5", "-r", "0.1", "-pct", "-s"], {}),
    ("flank", ["flank", "-i", D + "/flank/a.bed", "-g", D + "/flank/tiny.genome", "-b", "5"], {}),
    ("flank-l-s", ["flank", "-i", D + "/flank/a.bed", "-g", D + "/flank/tiny.genome", "-l", "3", "-r", "0", "-s"], {}),
    ("shift", ["shift", "-i", D + "/shift/a.bed", "-g", D + "/shift/tiny.genome", "-s", "5"], {}),
    ("shift-p-m-pct", ["shift", "-i", D + "/shift/a.bed", "-g", D + "/shift/tiny.genome", "-p", "0.5", "-m", "-0.2",
                       "-pct"], {}),
    ("makewindows-b", ["makewindows", "-b", D + "/makewindows/input.bed", "-w", "5000", "-s", "2000", "-i", "winnum"], {}),
    ("makewindows-g", ["makewindows", "-g", TINY, "-n", "7"], {}),
    ("makewindows-reverse", ["makewindows", "-b", D + "/makewindows/input.bed", "-n", "3", "-reverse", "-i", "srcwinnum"], {}),
    ("random", ["random", "-l", "100", "-n", "20", "-seed", "42", "-g", G19], {}),
    ("shuffle", ["shuffle", "-i", I + "a.bed", "-g", G19, "-seed", "7"], {}),
    ("shuffle-chrom-incl-excl", ["shuffle", "-i", D + "/shuffle/simrep.bed", "-g", G19, "-seed", "3", "-chrom",
                                 "-incl", D + "/shuffle/incl.bed"], {}),
    ("shuffle-excl", ["shuffle", "-i", I + "a.bed", "-g", G19, "-seed", "11", "-excl", D + "/shuffle/excl.bed",
                      "-noOverlapping"], {}),
    ("sample", ["sample", "-i", D + "/split/randData.bed", "-n", "10", "-seed", "5"], {}),
    ("sample-s", ["sample", "-i", I + "a.bed", "-n", "1", "-seed", "3", "-s", "forward"], {}),
    # statistics
    ("jaccard", ["jaccard", "-a", D + "/jaccard/a.bed", "-b", D + "/jaccard/b.bed"], {}),
    ("jaccard-split-bam", ["jaccard", "-a", D + "/jaccard/a.bam", "-b", D + "/jaccard/three_blocks_match.bam", "-bed"], {}),
    ("fisher", ["fisher", "-a", D + "/fisher/a.bed", "-b", D + "/fisher/b.bed", "-g", D + "/fisher/t.500.genome"], {}),
    ("fisher-m", ["fisher", "-a", D + "/fisher/a_merge.bed", "-b", D + "/fisher/b.bed", "-g", D + "/fisher/t.60.genome",
                  "-m"], {}),
    ("fisher-gff", ["fisher", "-b", D + "/fisher/test.bed", "-a", D + "/fisher/tumor.gff", "-g", D + "/fisher/dm6.fai"], {}),
    ("reldist", ["reldist", "-a", D + "/reldist/issue_711.a.bed", "-b", D + "/reldist/issue_711.b.bed"], {}),
    ("reldist-detail", ["reldist", "-a", D + "/reldist/issue_711.a.bed", "-b", D + "/reldist/issue_711.b.bed",
                        "-detail"], {}),
    ("summary", ["summary", "-i", I + "a.bed", "-g", G19], {}),
    ("annotate", ["annotate", "-i", I + "a.bed", "-files", I + "b.bed", I + "c.bed"], {}),
    ("annotate-counts-names", ["annotate", "-i", I + "a.bed", "-files", I + "b.bed", I + "c.bed", "-names", "B", "C",
                               "-counts", "-s"], {}),
    ("overlap", ["overlap", "-i", "ovl.txt", "-cols", "2,3,5,6"], {"setup": {"ovl.txt": OVL}}),
    # multiple files
    ("multiinter", ["multiinter", "-i", I + "a.bed", I + "b.bed", I + "c.bed"], {}),
    ("multiinter-header", ["multiinter", "-header", "-i", I + "a.bed", I + "b.bed", I + "c.bed"], {}),
    ("multiinter-header-names", ["multiinter", "-header", "-names", "A", "B", "C", "-i", I + "a.bed", I + "b.bed",
                                 I + "c.bed"], {}),
    ("multiinter-header-relative", ["multiinter", "-header", "-i", "a.bed", "sub/b.bed"],
     {"setup": {"a.bed": I + "a.bed", "sub/b.bed": I + "b.bed"}}),
    ("multiinter-header-native", ["multiinter", "-header", "-i", I + "a.bed", SP + "/b.bed"],
     {"setup": {"a b/b.bed": I + "b.bed"}, "native": True}),
    ("multiinter-empty-cluster", ["multiinter", "-header", "-empty", "-g", TINY, "-i", I + "a.bed", I + "b.bed",
                                  "-cluster"], {}),
    ("multiinter-nonascii", ["multiinter", "-header", "-i", NZ + "/\u0105.bed", NZ + "/b.bed"],
     {"setup": {"za\u017c\u00f3\u0142\u0107/\u0105.bed": I + "a.bed", "za\u017c\u00f3\u0142\u0107/b.bed": I + "b.bed"}}),
    ("unionbedg", ["unionbedg", "-i", "1.bg", "2.bg", "3.bg"],
     {"setup": {"1.bg": BG[0], "2.bg": BG[1], "3.bg": BG[2]}}),
    ("unionbedg-header", ["unionbedg", "-header", "-i", C + "/1.bg", C + "/2.bg", C + "/3.bg"],
     {"setup": {"1.bg": BG[0], "2.bg": BG[1], "3.bg": BG[2]}}),
    ("unionbedg-header-names-empty", ["unionbedg", "-header", "-names", "x", "y", "z", "-empty", "-g", "g.txt",
                                      "-filler", "N/A", "-i", "1.bg", "2.bg", "3.bg"],
     {"setup": {"1.bg": BG[0], "2.bg": BG[1], "3.bg": BG[2], "g.txt": GENOME2}}),
    ("unionbedg-one-file", ["unionbedg", "-i", "1.bg"], {"setup": {"1.bg": BG[0]}}),
    ("multicov", ["multicov", "-bams", I + "two_blocks_partial.bam", "-bed", I + "two_blocks_partial.bed"], {}),
    ("multicov-split", ["multicov", "-split", "-bams", I + "two_blocks_partial.bam", "-bed",
                        I + "two_blocks_partial.bed"], {}),
    # FASTA
    ("getfasta", ["getfasta", "-fi", FA, "-bed", "r.bed"], {"setup": R}),
    ("getfasta-fo", ["getfasta", "-fi", FA, "-bed", "r.bed", "-fo", "out.fa"], {"setup": R}),
    ("getfasta-s-name", ["getfasta", "-fi", FA, "-bed", D + "/getfasta/blocks.bed", "-s", "-name"], {}),
    ("getfasta-tab-split", ["getfasta", "-fi", FA, "-bed", D + "/getfasta/blocks.bed", "-split", "-tab", "-s"], {}),
    ("getfasta-bgzf-fasta", ["getfasta", "-fi", "t.fa.gz", "-bed", "r.bed"],
     {"setup": {"r.bed": REGIONS, "t.fa.gz": D + "/getfasta/t.fa.gz", "t.fa.gz.fai": D + "/getfasta/t.fa.gz.fai",
                "t.fa.gz.gzi": D + "/getfasta/t.fa.gz.gzi"}}),
    ("getfasta-makes-fai", ["getfasta", "-fi", "t.fa", "-bed", "r.bed"], {"setup": {"t.fa": FA, "r.bed": REGIONS}}),
    ("getfasta-crlf-fasta", ["getfasta", "-fi", "crlf.fa", "-bed", "r.bed"],
     {"setup": {"crlf.fa": crlf(FA), "r.bed": REGIONS}}),
    ("getfasta-rna-gff", ["getfasta", "-fi", "rna.fasta", "-bed", D + "/getfasta/rna.bed", "-rna"],
     {"setup": {"rna.fasta": D + "/getfasta/rna.fasta"}}),
    ("getfasta-missing-fasta", ["getfasta", "-fi", "nosuch.fa", "-bed", "r.bed"], {"setup": R}),
    ("getfasta-bed-stdin", ["getfasta", "-fi", FA, "-bed", "-", "-s", "-nameOnly"], {"stdin": REGIONS}),
    ("maskfasta-soft", ["maskfasta", "-fi", FA, "-bed", "r.bed", "-fo", "masked.fa", "-soft"], {"setup": R}),
    ("maskfasta", ["maskfasta", "-fi", FA, "-bed", "r.bed", "-fo", "masked.fa"], {"setup": R}),
    ("maskfasta-mc", ["maskfasta", "-fi", FA, "-bed", "r.bed", "-fo", "masked.fa", "-mc", "X"], {"setup": R}),
    ("nuc", ["nuc", "-fi", FA, "-bed", "r.bed"], {"setup": R}),
    ("nuc-s-seq-pattern", ["nuc", "-fi", FA, "-bed", D + "/getfasta/blocks.bed", "-s", "-seq", "-pattern", "AC", "-C"], {}),
    # BAM conversions
    ("bamtobed", ["bamtobed", "-i", I + "a.bam"], {}),
    ("bamtobed-bed12", ["bamtobed", "-bed12", "-i", I + "split.issue750.bam"], {}),
    ("bamtobed-split", ["bamtobed", "-split", "-i", I + "split.issue750.bam"], {}),
    ("bamtobed-bedpe", ["bamtobed", "-bedpe", "-i", I + "a_with_bothUnmapped.bam"], {}),
    ("bamtobed-cigar", ["bamtobed", "-cigar", "-i", I + "x.bam"], {}),
    ("bamtobed-cram", ["bamtobed", "-i", I + "a.cram"], {"env": {"CRAM_REFERENCE": I + "test_ref.fa"}}),
    ("bamtobed-stdin", ["bamtobed", "-i", "-"], {"stdin": "intersect/a.bam"}),
    ("bamtobed-stdin-empty", ["bamtobed", "-i", "-"], {"stdin": b""}),
    ("bamtobed-missing", ["bamtobed", "-i", "nosuch.bam"], {}),
    ("bamtofastq", ["bamtofastq", "-i", I + "split.issue750.bam", "-fq", "out.fq"], {}),
    ("bamtofastq-tags", ["bamtofastq", "-i", I + "a.bam", "-fq", "out.fq", "-tags"], {}),
    ("bedtobam", ["bedtobam", "-i", I + "a.bed", "-g", G19], {}),
    ("bedtobam-ubam-bed12", ["bedtobam", "-i", I + "blocks.bed12", "-g", G19, "-bed12", "-ubam"], {}),
    ("bedtobam-stdin", ["bedtobam", "-i", "-", "-g", G19], {"stdin": "intersect/a.bed"}),
    ("bedpetobam", ["bedpetobam", "-i", "p.bedpe", "-g", "g.txt"], {"setup": {"p.bedpe": BEDPE, "g.txt": GENOME2}}),
    ("bed12tobed6", ["bed12tobed6", "-i", D + "/bed12tobed6/three_blocks.bed"], {}),
    ("bed12tobed6-n", ["bed12tobed6", "-i", D + "/bed12tobed6/two_blocks.bed", "-n"], {}),
    ("tag", ["tag", "-i", I + "a.bam", "-files", I + "b.bed", "-labels", "B"], {}),
    ("tag-names-intervals", ["tag", "-i", I + "a.bam", "-files", I + "b.bed", "-names", "-tag", "YB"], {}),
    # pairs
    ("pairtobed", ["pairtobed", "-a", "p.bedpe", "-b", I + "a.bed"], {"setup": {"p.bedpe": BEDPE}}),
    ("pairtobed-type-neither", ["pairtobed", "-a", "p.bedpe", "-b", I + "a.bed", "-type", "neither"],
     {"setup": {"p.bedpe": BEDPE}}),
    ("pairtobed-abam", ["pairtobed", "-abam", I + "a_with_bothUnmapped.bam", "-b", I + "a.bed", "-bedpe",
                        "-type", "either"], {}),
    ("pairtopair", ["pairtopair", "-a", "p.bedpe", "-b", "q.bedpe"], {"setup": {"p.bedpe": BEDPE, "q.bedpe": BEDPE_B}}),
    ("pairtopair-type-either", ["pairtopair", "-a", "p.bedpe", "-b", "q.bedpe", "-type", "either", "-slop", "50"],
     {"setup": {"p.bedpe": BEDPE, "q.bedpe": BEDPE_B}}),
    # text tools
    ("groupby", ["groupby", "-i", D + "/groupby/values3.header.bed", "-c", "5", "-header"], {}),
    ("groupby-ops", ["groupby", "-i", D + "/groupby/test.bed", "-g", "1", "-c", "5,5,5,4", "-o",
                     "mean,stdev,median,collapse"], {}),
    ("groupby-full-ignorecase", ["groupby", "-i", D + "/groupby/values3_case.header.bed", "-c", "5", "-ignorecase",
                                 "-full"], {}),
    ("groupby-bam", ["groupby", "-i", D + "/groupby/gdc.bam", "-g", "1,3", "-c", "4", "-o", "mean"], {}),
    ("groupby-stdin", ["groupby", "-i", "-", "-g", "1", "-c", "2,3", "-o", "distinct,min"],
     {"stdin": "groupby/bug569_problem.txt"}),
    ("expand", ["expand", "-i", D + "/expand/expand.txt", "-c", "5,4"], {}),
    ("split-size", ["split", "-i", D + "/split/randData.bed", "-n", "3", "-p", "part", "-a", "size"], {}),
    ("split-simple", ["split", "-i", D + "/split/randData.bed", "-n", "4", "-p", C + "/sub dir/p", "-a", "simple"],
     {"setup": {"sub dir/.keep": b""}}),
    ("igv", ["igv", "-i", I + "a.bed"], {}),
    ("igv-path-name", ["igv", "-i", I + "a.bed", "-path", C + "/shots", "-sess", "s.xml", "-name", "-slop", "10",
                       "-img", "svg"], {}),
    ("links", ["links", "-i", I + "a.bed"], {}),
    ("links-native", ["links", "-i", I + "a.bed"], {"native": True}),
    ("links-base-org", ["links", "-i", I + "a.bed", "-base", "http://example.org", "-org", "mouse", "-db", "mm10"], {}),
]


# ---------------------------------------------------------------- dump
def path_forms(path):
    """Spellings of a directory in output: as given and resolved, both
    separators, and the drive forms C:\\x, C:/x, c:/x, /c/x."""
    forms = set()
    for p in {path, os.path.abspath(path), os.path.realpath(path)}:
        for q in {p, p.replace("\\", "/"), p.replace("/", "\\")}:
            forms.add(q)
            m = re.match(r"^([A-Za-z]):[\\/](.*)$", q)
            if m:
                d, rest = m.group(1), m.group(2)
                for dd in {d.upper(), d.lower()}:
                    forms |= {dd + ":/" + rest.replace("\\", "/"), dd + ":\\" + rest.replace("/", "\\")}
                forms.add("/" + d.lower() + "/" + rest.replace("\\", "/"))
    return forms


def norm(data, cwd, data_dir):
    """Bytes with the directories replaced by <CWD>/<DATA> (UTF-8, and on
    Windows also the ANSI code page), decoded as latin-1 (every byte kept)."""
    pairs = []
    for path, tag in ((cwd, b"<CWD>"), (data_dir, b"<DATA>")):
        for f in path_forms(path):
            for enc in ("utf-8", "mbcs") if os.name == "nt" else ("utf-8",):
                try:
                    pairs.append((f.encode(enc, "strict"), tag))
                except (UnicodeEncodeError, LookupError):
                    pass
    for form, tag in sorted(set(pairs), key=lambda ft: -len(ft[0])):
        data = data.replace(form, tag)
    return data.decode("latin-1")


def bam_text(raw):
    """Readable, length-independent form of uncompressed BAM: header text,
    references, then one hex line per alignment record."""
    i = 4
    l_text = int.from_bytes(raw[i:i + 4], "little"); i += 4
    lines = [b"#BAM header:\n" + raw[i:i + l_text].rstrip(b"\0")]; i += l_text
    n_ref = int.from_bytes(raw[i:i + 4], "little"); i += 4
    for _ in range(n_ref):
        l_name = int.from_bytes(raw[i:i + 4], "little"); i += 4
        name = raw[i:i + l_name].rstrip(b"\0"); i += l_name
        lines.append(b"#ref " + name + b" " + str(int.from_bytes(raw[i:i + 4], "little")).encode()); i += 4
    while i + 4 <= len(raw):
        size = int.from_bytes(raw[i:i + 4], "little")
        lines.append(b"#rec " + raw[i + 4:i + 4 + size].hex().encode()); i += 4 + size
    if i != len(raw):
        lines.append(b"#trailing " + raw[i:].hex().encode())
    return b"\n".join(lines) + b"\n"


def content(raw):
    """(kind, decoded bytes): gzip/BGZF decompressed, BAM made readable."""
    kind = "plain"
    if raw[:2] == b"\x1f\x8b":
        kind = "bgzf" if raw[3:4] == b"\x04" and raw[12:14] == b"BC" else "gzip"
        try:
            raw = gzip.decompress(raw)
        except (OSError, EOFError) as e:
            return kind + "+corrupt(" + str(e) + ")", raw
    if raw[:4] == b"CRAM":              # blocks inside are deflated; no decoder here
        return kind + "+cram", raw
    if raw[:4] == b"BAM\1":
        kind += "+bam"
        try:
            raw = bam_text(raw)
        except Exception as e:          # keep the bytes, say why
            kind += "+unparsed(" + str(e) + ")"
    return kind, raw


def make_inputs(spec, cwd):
    for dst, src in spec.items():
        dst = os.path.join(cwd, *dst.split("/"))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if isinstance(src, bytes):
            data = src
        elif isinstance(src, tuple):
            data = open(subst(src[1], cwd), "rb").read()
            data = data.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n") if src[0] == "crlf" else gzip.compress(data, mtime=0)
        else:
            data = open(subst(src, cwd), "rb").read()
        with open(dst, "wb") as f:
            f.write(data)


def subst(arg, cwd):
    return arg.replace(D, DATA).replace(C, cwd)


def snapshot(cwd):
    files = {}
    for root, _, fs in os.walk(cwd):
        for f in fs:
            full = os.path.join(root, f)
            files[os.path.relpath(full, cwd).replace(os.sep, "/")] = open(full, "rb").read()
    return files


def run(case, cwd):
    name, args, opt = case
    shutil.rmtree(cwd, ignore_errors=True)
    os.makedirs(cwd)
    make_inputs(opt.get("setup", {}), cwd)
    before = snapshot(cwd)
    env = dict(os.environ)
    for k, v in opt.get("env", {}).items():
        env[k] = subst(v, cwd)
    stdin = opt.get("stdin")
    if stdin is None:
        inp = subprocess.DEVNULL
    elif isinstance(stdin, bytes):
        inp = stdin
    else:
        inp = open(os.path.join(DATA, *stdin.split("/")), "rb").read()
    kw = {"input": inp} if isinstance(inp, bytes) else {"stdin": inp}
    argv = [subst(a, cwd) for a in args]
    if opt.get("native"):
        # MSYS2's Python may use / as os.sep, so spell the backslashes out
        argv = [os.path.normpath(a).replace("/", "\\") if os.name == "nt" and a != s else a
                for a, s in zip(argv, args)]
    p = subprocess.run([EXE] + argv, cwd=cwd, capture_output=True, env=env, **kw)
    files = {}
    for rel, data in sorted(snapshot(cwd).items()):
        if before.get(rel) != data:
            kind, data = content(data)
            files[rel] = [kind, norm(data, cwd, DATA)]
    kind, out = content(p.stdout)
    return {"exit": p.returncode, "stdout": [kind, norm(out, cwd, DATA)],
            "stderr": norm(p.stderr, cwd, DATA), "files": files}


# ---------------------------------------------------------------- cmp
# Text that differs by design and is replaced before comparing (both sides):
# the version (the macOS build is a git checkout, `git describe` adds
# -N-gHASH; a build from an exported tree prints the plain release), and the
# separator inside a path after <CWD>/<DATA> (Windows joins with \ itself).
XPLAT_TEXT = [
    (r"v2\.\d+\.\d+(?:-\d+-g[0-9a-f]+)?(?:-dirty)?", "v<VERSION>"),
]


def xnorm(text):
    for pattern, repl in XPLAT_TEXT:
        text = re.sub(pattern, repl, text)
    return re.sub(r"<(CWD|DATA)>[^\t\r\n\"'<>()]*", lambda m: m.group(0).replace("\\", "/"), text)


# Known differences between the Windows and the macOS build: case -> reason.
EXPECTED_XPLAT = {
    "intersect-nonascii-path": "Windows passes arguments in the ANSI code page unless the .exe embeds a UTF-8"
                               " manifest (activeCodePage); a plain build has none, the ngs-tools payload does",
    "multiinter-nonascii": "same as intersect-nonascii-path",
    "intersect-cram": "CRAM output: its blocks are deflated by zlib (macOS) or libdeflate (Windows) and"
                      " cannot be decoded here; intersect-cram-bed checks the same records as text",
    "split-size": "upstream std::sort of equal-sized records is unstable; libstdc++ (GCC, UCRT64 and Linux)"
                  " orders ties unlike libc++ (clang: macOS, CLANGARM64), so records land in other files",
}


def show(text):
    """Lines with CR and other control bytes made visible."""
    return [repr(l)[1:-1] for l in text.split("\n")]


def diff(a, b, label):
    d = list(difflib.unified_diff(show(a), show(b), "A", "B", lineterm="", n=1))
    return [f"    {label}:"] + ["      " + x for x in d[:20]] + (["      ..."] if len(d) > 20 else [])


def compare(a, b):
    problems = []
    if a["exit"] != b["exit"]:
        problems.append(f"    exit: {a['exit']} vs {b['exit']}")
    streams = [("stdout", a["stdout"], b["stdout"]), ("stderr", ["plain", a["stderr"]], ["plain", b["stderr"]])]
    streams += [(f"file {f}", a["files"].get(f), b["files"].get(f)) for f in sorted(set(a["files"]) | set(b["files"]))]
    for label, x, y in streams:
        if x is None or y is None:
            problems.append(f"    {label}: only in {'B' if x is None else 'A'}")
            continue
        if x[0] != y[0]:
            problems.append(f"    {label}: kind {x[0]} vs {y[0]}")
        tx, ty = xnorm(x[1]), xnorm(y[1])
        if tx != ty and "cram" in x[0] + y[0]:
            problems.append(f"    {label}: CRAM bytes differ ({len(tx)} vs {len(ty)} bytes; not decoded)")
        elif tx != ty:
            problems += diff(tx, ty, label)
    return problems


def main():
    global EXE, DATA
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "dump" and len(sys.argv) in (6, 7):
        EXE, DATA, work = (os.path.abspath(p) for p in sys.argv[2:5])
        only = sys.argv[6] if len(sys.argv) > 6 else ""
        out = {}
        for case in CASES:
            if only in case[0]:
                out[case[0]] = run(case, os.path.join(work, case[0]))
                print("ran", case[0], "exit", out[case[0]]["exit"], flush=True)
        with open(sys.argv[5], "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1, ensure_ascii=True)
        print(len(out), "case(s) written to", sys.argv[5])
        return
    if mode != "cmp" or len(sys.argv) != 4:
        sys.exit(__doc__)
    a, b = (json.load(open(p, encoding="utf-8")) for p in sys.argv[2:4])
    counts = {"SAME": 0, "EXPECTED-DIFF": 0, "DIFF": 0}
    for name in a:
        if name not in b:
            print(f"MISSING  {name} (only in A)")
            continue
        problems = compare(a[name], b[name])
        status = "SAME" if not problems else ("EXPECTED-DIFF" if name in EXPECTED_XPLAT else "DIFF")
        counts[status] += 1
        print(f"{status}  {name}" + (f"  ({EXPECTED_XPLAT[name]})" if status == "EXPECTED-DIFF" else ""))
        if status == "DIFF":
            print("\n".join(problems))
    for name in b:
        if name not in a:
            print(f"MISSING  {name} (only in B)")
    print(f"\n{counts['SAME']} same, {counts['EXPECTED-DIFF']} expected, {counts['DIFF']} different")
    sys.exit(1 if counts["DIFF"] else 0)


if __name__ == "__main__":
    main()
