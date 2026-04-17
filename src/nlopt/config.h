/* config.h — Portable configuration for NLopt bundled in SMAD R package.
 *
 * This replaces the autoconf-generated config.h with compile-time platform
 * detection so that the same source tree builds on:
 *   - Linux   x86_64 / arm64   (GCC / Clang)
 *   - macOS   x86_64 / arm64   (Apple Clang)
 *   - Windows x86_64           (MinGW-w64 via Rtools)
 *
 * Most feature macros are also passed via PKG_CFLAGS in Makevars /
 * Makevars.win, but this header provides fallbacks and platform-specific
 * type-size detection that cannot be done from Makevars alone.
 */

#ifndef NLOPT_CONFIG_H
#define NLOPT_CONFIG_H

/* ===== Version ========================================================== */
#define MAJOR_VERSION 2
#define MINOR_VERSION 3
#define BUGFIX_VERSION 0
#define PACKAGE "nlopt"
#define PACKAGE_BUGREPORT "stevenj@alum.mit.edu"
#define PACKAGE_NAME "nlopt"
#define PACKAGE_STRING "nlopt 2.3"
#define PACKAGE_TARNAME "nlopt"
#define PACKAGE_URL ""
#define PACKAGE_VERSION "2.3"
#define VERSION "2.3"

/* ===== Standard C headers (always available) ============================= */
#define STDC_HEADERS 1

#ifndef HAVE_STDINT_H
#define HAVE_STDINT_H 1
#endif

#define HAVE_STDLIB_H 1
#define HAVE_STRING_H 1
#define HAVE_STRINGS_H 1
#define HAVE_MEMORY_H 1

/* ===== Math functions (universally available) ============================ */
#ifndef HAVE_COPYSIGN
#define HAVE_COPYSIGN 1
#endif

#ifndef HAVE_ISINF
#define HAVE_ISINF 1
#endif

#ifndef HAVE_ISNAN
#define HAVE_ISNAN 1
#endif

#define HAVE_LIBM 1

/* ===== Time ============================================================= */
#ifndef HAVE_TIME
#define HAVE_TIME 1
#endif

/* ===== Platform-specific: POSIX vs Windows ============================== */
#ifdef _WIN32
  /* ---- Windows (MinGW-w64 / MSVC) ------------------------------------- */
  /* No unistd.h, no dlfcn.h, no gettimeofday, no getopt.h */
  /* getpid is available as _getpid */
  #define HAVE_GETPID 1

  /* sizeof(unsigned long) is 4 on Windows (LLP64 model) */
  #define SIZEOF_UNSIGNED_INT  4
  #define SIZEOF_UNSIGNED_LONG 4
  #define HAVE_UINT32_T 1

  /* Thread-local storage */
  #if defined(__GNUC__)
    /* MinGW GCC supports __thread */
    #define THREADLOCAL __thread
  #elif defined(_MSC_VER)
    #define THREADLOCAL __declspec(thread)
  #else
    #define THREADLOCAL
  #endif

#else
  /* ---- POSIX: Linux / macOS (any architecture) ------------------------- */
  #ifndef HAVE_UNISTD_H
  #define HAVE_UNISTD_H 1
  #endif

  #define HAVE_DLFCN_H 1
  #define HAVE_GETOPT_H 1
  #define HAVE_GETPID 1
  #define HAVE_SYS_STAT_H 1
  #define HAVE_SYS_TYPES_H 1

  #ifndef HAVE_GETTIMEOFDAY
  #define HAVE_GETTIMEOFDAY 1
  #endif

  #define TIME_WITH_SYS_TIME 1

  /* sizeof(unsigned long) is 8 on LP64 (Linux/macOS 64-bit) */
  #define SIZEOF_UNSIGNED_INT  4
  #define SIZEOF_UNSIGNED_LONG 8
  #define HAVE_UINT32_T 1

  /* qsort_r: available on glibc (Linux) and macOS */
  #if defined(__GLIBC__) || defined(__APPLE__)
    #define HAVE_QSORT_R 1
  #endif

  /* gettid syscall: Linux only */
  #ifdef __linux__
    #define HAVE_GETTID_SYSCALL 1
  #endif

  /* Thread-local storage: GCC/Clang */
  #define THREADLOCAL __thread

#endif /* _WIN32 */

/* ===== Optional features (disabled) ===================================== */
/* #undef DEBUG */
/* #undef WITH_CXX */
/* #undef WITH_NOCEDAL */

/* ===== Inline =========================================================== */
#ifndef __cplusplus
/* #undef inline */
#endif

/* ===== R compatibility macros =========================================== */
#include <R.h>
#include <Rinternals.h>
#define printf Rprintf
#define R_STDOUT ((FILE*)1)
#define R_STDERR ((FILE*)2)
#undef stdout
#define stdout R_STDOUT
#undef stderr
#define stderr R_STDERR
#define fprintf(f, ...) ( ((f) == R_STDERR) ? REprintf(__VA_ARGS__) : ( ((f) == R_STDOUT) ? Rprintf(__VA_ARGS__) : (f) == NULL ? Rprintf(__VA_ARGS__) : fprintf(f, __VA_ARGS__) ) )
#define putchar Rputchar
#define puts(s) Rprintf("%s\n", (s))
#ifndef exit
#define exit(x) error("Exit called with code %d", (x))
#endif

/* Suppress deprecation warnings in NLopt sources */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#endif /* NLOPT_CONFIG_H */
