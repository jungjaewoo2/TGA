/* include/config.h.  Generated from config.h.in by configure.  */
/* include/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to recognise all pointers to the interior of objects. */
#define ALL_INTERIOR_POINTERS 1

/* Define to enable atomic uncollectible allocation. */
#define ATOMIC_UNCOLLECTABLE 1

/* See doc/README.macros. */
/* #undef DARWIN_DONT_PARSE_STACK */

/* Define to force debug headers on all objects. */
/* #undef DBG_HDRS_ALL */

/* Define to enable support for DB/UX threads. */
/* #undef DGUX_THREADS */

/* Define to enable eCos target support. */
/* #undef ECOS */

/* Wine getenv may not return NULL for missing entry. */
/* #undef EMPTY_GETENV_RESULTS */

/* Define to enable alternative finalization interface. */
#define ENABLE_DISCLAIM 1

/* Define to support IBM AIX threads. */
/* #undef GC_AIX_THREADS */

/* Define to enable internal debug assertions. */
/* #undef GC_ASSERTIONS */

/* Define to support Darwin pthreads. */
/* #undef GC_DARWIN_THREADS */

/* Define to enable support for DB/UX threads on i386. */
/* #undef GC_DGUX386_THREADS */

/* Define to build dynamic libraries with only API symbols exposed. */
/* #undef GC_DLL */

/* Define to support FreeBSD pthreads. */
/* #undef GC_FREEBSD_THREADS */

/* Define to include support for gcj. */
#define GC_GCJ_SUPPORT 1

/* Define to support GNU pthreads. */
/* #undef GC_GNU_THREADS */

/* Define if backtrace information is supported. */
/* #undef GC_HAVE_BUILTIN_BACKTRACE */

/* Define to support HP/UX 11 pthreads. */
/* #undef GC_HPUX_THREADS */

/* Enable Win32 DllMain-based approach of threads registering. */
/* #undef GC_INSIDE_DLL */

/* Define to support Irix pthreads. */
/* #undef GC_IRIX_THREADS */

/* Define to support pthreads on Linux. */
/* #undef GC_LINUX_THREADS */

/* Define to support NetBSD pthreads. */
/* #undef GC_NETBSD_THREADS */

/* Define to support OpenBSD pthreads. */
/* #undef GC_OPENBSD_THREADS */

/* Define to support Tru64 pthreads. */
/* #undef GC_OSF1_THREADS */

/* Read environment variables from the GC 'env' file. */
/* #undef GC_READ_ENV_FILE */

/* Define to support rtems-pthreads. */
/* #undef GC_RTEMS_PTHREADS */

/* Define to support Solaris pthreads. */
/* #undef GC_SOLARIS_THREADS */

/* Define to support platform-specific threads. */
/* #undef GC_THREADS */

/* Explicitly prefix exported/imported WINAPI symbols with '_'. */
/* #undef GC_UNDERSCORE_STDCALL */

/* Force the GC to use signals based on SIGRTMIN+k. */
/* #undef GC_USESIGRT_SIGNALS */

/* See doc/README.macros. */
/* #undef GC_USE_DLOPEN_WRAP */

/* The major version number of this GC release. */
#define GC_VERSION_MAJOR 7

/* The micro version number of this GC release. */
#define GC_VERSION_MICRO 4

/* The minor version number of this GC release. */
#define GC_VERSION_MINOR 4

/* Define to support pthreads-win32. */
/* #undef GC_WIN32_PTHREADS */

/* Define to support Win32 threads. */
/* #undef GC_WIN32_THREADS */

/* Define to install pthread_atfork() handlers by default. */
/* #undef HANDLE_FORK */

/* Define to use 'dladdr' function. */
#define HAVE_DLADDR 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* See doc/README.macros. */
#define JAVA_FINALIZATION 1

/* Define to save back-pointers in debugging headers. */
/* #undef KEEP_BACK_PTRS */

/* Define to optimize for large heaps or root sets. */
/* #undef LARGE_CONFIG */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* See doc/README.macros. */
/* #undef MAKE_BACK_GRAPH */

/* Number of GC cycles to wait before unmapping an unused block. */
/* #undef MUNMAP_THRESHOLD */

/* Define to not use system clock (cross compiling). */
/* #undef NO_CLOCK */

/* Disable debugging, like GC_dump and its callees. */
/* #undef NO_DEBUGGING */

/* Define to make the collector not allocate executable memory by default. */
#define NO_EXECUTE_PERMISSION 1

/* Prohibit installation of pthread_atfork() handlers. */
/* #undef NO_HANDLE_FORK */

/* Name of package */
#define PACKAGE "gc"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "bdwgc@lists.opendylan.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "gc"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "gc 7.4.4"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "gc"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "7.4.4"

/* Define to enable parallel marking. */
/* #undef PARALLEL_MARK */

/* If defined, redirect free to this function. */
/* #undef REDIRECT_FREE */

/* If defined, redirect malloc to this function. */
/* #undef REDIRECT_MALLOC */

/* If defined, redirect GC_realloc to this function. */
/* #undef REDIRECT_REALLOC */

/* The number of caller frames saved when allocating with the debugging API.
   */
/* #undef SAVE_CALL_COUNT */

/* Shorten the headers to minimize object size at the expense of checking for
   writes past the end (see doc/README.macros). */
/* #undef SHORT_DBG_HDRS */

/* Define to tune the collector for small heap sizes. */
/* #undef SMALL_CONFIG */

/* See the comment in gcconfig.h. */
/* #undef SOLARIS25_PROC_VDB_BUG_FIXED */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to work around a Solaris 5.3 bug (see dyn_load.c). */
/* #undef SUNOS53_SHARED_LIB */

/* Define to enable thread-local allocation optimization. */
/* #undef THREAD_LOCAL_ALLOC */

/* Use Unicode (W) variant of Win32 API instead of ASCII (A) one. */
/* #undef UNICODE */

/* Define to use of compiler-support for thread-local variables. */
/* #undef USE_COMPILER_TLS */

/* Define to use mmap instead of sbrk to expand the heap. */
/* #undef USE_MMAP */

/* Define to return memory to OS with munmap calls (see doc/README.macros). */
/* #undef USE_MUNMAP */

/* Define to use Win32 VirtualAlloc (instead of sbrk or mmap) to expand the
   heap. */
/* #undef USE_WINALLOC */

/* Version number of package */
#define VERSION "7.4.4"

/* The POSIX feature macro. */
/* #undef _POSIX_C_SOURCE */

/* Indicates the use of pthreads (NetBSD). */
/* #undef _PTHREADS */

/* Required define if using POSIX threads. */
/* #undef _REENTRANT */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
