#pragma once

#define _FILE_OFFSET_BITS 64

#ifdef WIN32
#include <stdint.h> // bedtools-crossplatform: int64_t
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef int64_t off_type; // bedtools-crossplatform: __int64_t is not a MinGW type
#else
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off_t off_type;
#endif
