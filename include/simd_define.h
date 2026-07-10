#ifndef SIMD_DEFINE
#define SIMD_DEFINE

/*
 * SIMD type definitions and intrinsics compatibility layer.
 *
 * Three paths are supported:
 * 1. NBODY_USE_SIMDE: Use SIMDe (SIMD Everywhere) for ARM/cross-platform compatibility
 * 2. __USE_INTEL: Use Intel intrinsics directly (requires Intel compiler or compatible)
 * 3. __USE_GNU (default): Use GCC vector extensions (x86 GCC only)
 */

/* SIMDe path - for ARM and cross-platform compatibility */
#ifdef NBODY_USE_SIMDE

#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/x86/sse.h>
#include <simde/x86/sse2.h>
#include <simde/x86/sse3.h>
#include <simde/x86/sse4.1.h>

/*
 * Wrapper structs for SIMD types to enable operator overloading.
 * SIMDe types are unions, not classes, so we cannot overload operators directly.
 * These wrappers provide implicit conversion to/from the underlying SIMDe types
 * and support bitcast conversions between different SIMD types.
 */
struct v4sf;
struct v2df;
struct v4si;

struct v2df {
    __m128d v;
    v2df() : v(_mm_setzero_pd()) {}
    v2df(__m128d x) : v(x) {}
    operator __m128d() const { return v; }
    /* Bitcast from v4sf - reinterprets bits as double */
    explicit v2df(const v4sf& x);
};

struct v4sf {
    __m128 v;
    v4sf() : v(_mm_setzero_ps()) {}
    v4sf(__m128 x) : v(x) {}
    operator __m128() const { return v; }
    /* Bitcast from v2df - reinterprets bits as float */
    explicit v4sf(const v2df& x) : v(_mm_castpd_ps(x.v)) {}
};

/* Deferred implementation of v2df constructor from v4sf */
inline v2df::v2df(const v4sf& x) : v(_mm_castps_pd(x.v)) {}

struct v4si {
    __m128i v;
    v4si() : v(_mm_setzero_si128()) {}
    v4si(__m128i x) : v(x) {}
    operator __m128i() const { return v; }
};

/* Note: AVX types (v4df, v8sf) are not defined in SIMDe SSE mode.
   AVX is not supported on ARM with SIMDe - only SSE mode is available.
   Code using AVX types will need to be compiled separately for x86. */

/* SSE operator overloads */
inline v2df operator + (const v2df& a, const v2df& b) { return _mm_add_pd(a,b); }
inline v2df operator - (const v2df& a, const v2df& b) { return _mm_sub_pd(a,b); }
inline v2df operator * (const v2df& a, const v2df& b) { return _mm_mul_pd(a,b); }
inline v2df operator / (const v2df& a, const v2df& b) { return _mm_div_pd(a,b); }
inline v4sf operator + (const v4sf& a, const v4sf& b) { return _mm_add_ps(a,b); }
inline v4sf operator - (const v4sf& a, const v4sf& b) { return _mm_sub_ps(a,b); }
inline v4sf operator * (const v4sf& a, const v4sf& b) { return _mm_mul_ps(a,b); }
inline v4sf operator / (const v4sf& a, const v4sf& b) { return _mm_div_ps(a,b); }
inline v4si operator - (const v4si& a, const v4si& b) { return _mm_sub_epi32(a,b); }

/* Compound assignment operators */
inline v4sf& operator += (v4sf& a, const v4sf& b) { a = a + b; return a; }
inline v4sf& operator -= (v4sf& a, const v4sf& b) { a = a - b; return a; }
inline v4sf& operator *= (v4sf& a, const v4sf& b) { a = a * b; return a; }
inline v2df& operator += (v2df& a, const v2df& b) { a = a + b; return a; }
inline v2df& operator -= (v2df& a, const v2df& b) { a = a - b; return a; }
inline v2df& operator *= (v2df& a, const v2df& b) { a = a * b; return a; }

/* Prefetch - use GCC builtin or platform-specific */
#ifndef __builtin_prefetch
#define __builtin_prefetch(p,rw,i) ((void)0)
#endif

#define REP4(x) _mm_set1_ps(x)

/* Intel intrinsics path */
#elif defined(__USE_INTEL)

#include "immintrin.h"
#include "xmmintrin.h"
#include "emmintrin.h"
#include "pmmintrin.h"
#include "smmintrin.h"

typedef __m256d v4df;
typedef __m256  v8sf;
typedef __m128d v2df;
typedef __m128  v4sf;
typedef __m128i v4si;
// SSE
#define __builtin_prefetch(p,rw,i)                 _mm_prefetch(p,i)
#define __builtin_ia32_shufps(a,b,imm)             _mm_shuffle_ps(a,b,imm)
#define __builtin_ia32_cvtpd2ps(a)                 _mm_cvtpd_ps(a)
#define __builtin_ia32_minpd(a,b)                  _mm_min_pd(a,b)
#define __builtin_ia32_unpckhpd(a,b)               _mm_unpackhi_pd(a,b)
#define __builtin_ia32_unpcklpd(a,b)               _mm_unpacklo_pd(a,b)
#define __builtin_ia32_loaddqu(mem_addr)           _mm_loaddup_pd(mem_addr)
//#define __builtin_ia32_vec_ext_v2df(a,imm)         _mm_store_sd(a,imm)
// AVX
#define __builtin_ia32_rsqrtps256(a)               _mm256_rsqrt_ps(a)
#define __builtin_ia32_vinsertf128_ps256(a,b,imm)  _mm256_insertf128_ps(a,b,imm)
#define __builtin_ia32_cvtpd2ps256(a)              _mm256_cvtpd_ps(a)
#define __builtin_ia32_cvtps2pd256(a)              _mm256_cvtps_pd(a)
#define __builtin_ia32_movntps256(mem_addr,a)      _mm256_stream_ps(mem_addr,a)
#define __builtin_ia32_unpcklps256(a, b)           _mm256_unpacklo_ps(a,b)
#define __builtin_ia32_unpckhps256(a, b)           _mm256_unpackhi_ps(a,b)
#define __builtin_ia32_shufps256(a,b,imm)          _mm256_shuffle_ps(a,b,imm)
#define __builtin_ia32_vextractf128_ps256(a,imm)   _mm256_extractf128_ps(a,imm)
#define __builtin_ia32_haddpd256(a,b)              _mm256_hadd_pd(a,b)
#define __builtin_ia32_vextractf128_pd256(a,imm)   _mm256_extractf128_pd(a,imm)
#define __builtin_ia32_minpd256(a,b)               _mm256_min_pd(a,b)
//#define __builtin_ia32_movntdq(mem_addr,a)         _mm_stream_si128(mem_addr,a)

// SSE
inline v2df operator + (const v2df& a, const v2df& b) { return _mm_add_pd(a,b); }
inline v2df operator - (const v2df& a, const v2df& b) { return _mm_sub_pd(a,b); }
inline v2df operator * (const v2df& a, const v2df& b) { return _mm_mul_pd(a,b); }
inline v2df operator / (const v2df& a, const v2df& b) { return _mm_div_pd(a,b); }
inline v4sf operator + (const v4sf& a, const v4sf& b) { return _mm_add_ps(a,b); }
inline v4sf operator - (const v4sf& a, const v4sf& b) { return _mm_sub_ps(a,b); }
inline v4sf operator * (const v4sf& a, const v4sf& b) { return _mm_mul_ps(a,b); }
inline v4sf operator / (const v4sf& a, const v4sf& b) { return _mm_div_ps(a,b); }
inline v4si operator - (const v4si& a, const v4si& b) { return _mm_sub_epi32(a,b); }

// AVX
inline v4df operator + (const v4df& a, const v4df& b) { return _mm256_add_pd(a,b); }
inline v4df operator - (const v4df& a, const v4df& b) { return _mm256_sub_pd(a,b); }
inline v4df operator * (const v4df& a, const v4df& b) { return _mm256_mul_pd(a,b); }
inline v4df operator / (const v4df& a, const v4df& b) { return _mm256_div_pd(a,b); }
inline v4df& operator += (v4df &a, const v4df& b) { a = _mm256_add_pd(a,b); return a; }
inline v8sf operator + (const v8sf& a, const v8sf& b) { return _mm256_add_ps(a,b); }
inline v8sf operator - (const v8sf& a, const v8sf& b) { return _mm256_sub_ps(a,b); }
inline v8sf operator * (const v8sf& a, const v8sf& b) { return _mm256_mul_ps(a,b); }
inline v8sf operator / (const v8sf& a, const v8sf& b) { return _mm256_div_ps(a,b); }
inline v8sf& operator += (v8sf &a, const v8sf& b) { a = _mm256_add_ps(a,b); return a; }

#define REP4(x) {x,x,x,x}
#define REP8(x) {x,x,x,x,x,x,x,x}

/* GCC path - default for x86 GCC */
#else

#ifndef __USE_GNU
#define __USE_GNU
#endif

// SSE
typedef float  v4sf __attribute__((vector_size(16)));
typedef double v2df __attribute__((vector_size(16)));
//typedef int    v4si __attribute__((vector_size(16)));
//typedef long long v2di __attribute__ ((__vector_size__ (16)));
// AVX
typedef float  v8sf __attribute__((vector_size(32)));
typedef double v4df __attribute__((vector_size(32)));

#define REP4(x) {x,x,x,x}
#define REP8(x) {x,x,x,x,x,x,x,x}

#endif /* NBODY_USE_SIMDE / __USE_INTEL / __USE_GNU */

#endif /* SIMD_DEFINE */
