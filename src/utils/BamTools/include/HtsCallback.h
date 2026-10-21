// bedtools-crossplatform: the stream-callback I/O that bedtools used to add to
// its bundled htslib 1.9 (hfile.c/hts.c). It now lives here, built on htslib's
// exported backend API (hfile_init/hfile_destroy, hts_hopen), so bedtools can
// use an unmodified htslib.
#ifndef BEDTOOLS_HTS_CALLBACK_H
#define BEDTOOLS_HTS_CALLBACK_H

#include <sys/types.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#ifdef __cplusplus
extern "C" {
#endif

/// Defines the operations that used by an IO file
typedef struct hFILE_callback_ops {
	void* cb_data;
	ssize_t (*read)(void* cb_data, void* buf, size_t sz);
	ssize_t (*write)(void* cb_data, const void* buf, size_t sz);
	off_t (*seek)(void* cb_data, off_t ofs, int whence);
	int (*flush)(void* cb_data);
	int (*close)(void* cb_data);
} hFILE_callback_ops;

/// Open an hFILE whose reads/writes go through the callbacks in ops.
hFILE *hopen_callback(hFILE_callback_ops ops, const char* mode);

/// Open a SAM/BAM/CRAM/VCF/BCF/etc from stream callbacks (see hts_open for mode).
htsFile *hts_open_callback(const char* fn, hFILE_callback_ops* ops, const char* mode);

#ifdef __cplusplus
}
#endif

#endif
