// bedtools-crossplatform: stream-callback hFILE backend, moved out of the
// bundled htslib (see HtsCallback.h). Same behaviour as bedtools' former
// hfile.c additions, plus a NULL check after hfile_init.
#include <errno.h>
#include "HtsCallback.h"
#include "hfile_internal.h"

typedef struct {
	hFILE base;
	hFILE_callback_ops ops;
} hFILE_cb;

static ssize_t cb_read(hFILE *fpv, void *buffer, size_t nbytes)
{
	hFILE_cb* hcb = (hFILE_cb*)fpv;
	ssize_t ret;

	if(hcb->ops.read == NULL)
	{
		errno = EBADF;
		return -1;
	}

	do {
		ret = hcb->ops.read(hcb->ops.cb_data, buffer, nbytes);
	} while(ret < 0 && errno == EINTR);

	return ret;
}

static ssize_t cb_write(hFILE* fpv, const void* buffer, size_t nbytes)
{
	hFILE_cb* hcb = (hFILE_cb*)fpv;
	ssize_t ret;

	if(hcb->ops.write == NULL)
	{
		errno = EBADF;
		return -1;
	}

	do {
		ret = hcb->ops.write(hcb->ops.cb_data, buffer, nbytes);
	} while(ret < 0 && errno == EINTR);

	return ret;
}

static off_t cb_seek(hFILE *fpv, off_t offset, int whence)
{
	hFILE_cb* hcb = (hFILE_cb*)fpv;
	if(hcb->ops.seek)
		return hcb->ops.seek(hcb->ops.cb_data, offset, whence);

	errno = ESPIPE;
	return -1;
}

static int cb_flush(hFILE* fpv)
{
	hFILE_cb* hcb = (hFILE_cb*)fpv;

	return hcb->ops.flush ? hcb->ops.flush(hcb->ops.cb_data) : 0;
}

static int cb_close(hFILE* fpv)
{
	hFILE_cb* hcb = (hFILE_cb*)fpv;

	return hcb->ops.close ? hcb->ops.close(hcb->ops.cb_data) : 0;
}

static const struct hFILE_backend cb_backend =
{
	cb_read, cb_write, cb_seek, cb_flush, cb_close
};

extern "C" hFILE *hopen_callback(hFILE_callback_ops ops, const char* mode)
{
	hFILE_cb *ret = (hFILE_cb*) hfile_init(sizeof(*ret), mode, 0);
	if(ret == NULL) return NULL;

	ret->ops = ops;
	ret->base.backend = &cb_backend;

	return (hFILE*)ret;
}

extern "C" htsFile *hts_open_callback(const char* fn, hFILE_callback_ops* ops, const char* mode)
{
	if(NULL == ops) return NULL;
	hFILE* fp = hopen_callback(*ops, mode);
	if(NULL == fp) return NULL;
	htsFile* ret = hts_hopen(fp, fn ? fn : "-", mode);
	if(NULL == ret) hclose_abruptly(fp);
	return ret;
}
