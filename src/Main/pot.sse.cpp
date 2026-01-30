#include <cstdio>
#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
#include "simd_define.h"

/* Include SSE intrinsics if using SIMDe or native x86 */
#ifdef NBODY_USE_SIMDE
/* SIMDe headers already included via simd_define.h */
#elif !defined(__USE_INTEL)
#include <xmmintrin.h>  /* SSE */
#include <emmintrin.h>  /* SSE2 */
#include <pmmintrin.h>  /* SSE3 */
#include <smmintrin.h>  /* SSE4.1 */
#endif

static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

//typedef float v4sf __attribute__ ((vector_size(16)));
static inline v4sf v4sf_rsqrt(v4sf x){
	v4sf y = _mm_rsqrt_ps(x);
	return (_mm_set1_ps(-0.5f) * y) * 
			(x*y*y + _mm_set1_ps(-3.f));
}

/* Helper function to extract float from v4sf at specific index */
static inline float extract_ps(v4sf v, int index){
	/* Use shuffle to move desired element to position 0, then extract */
	union { v4sf v; float f[4]; } u;
	u.v = v;
	return u.f[index];
}

static inline void pot_reduce(v4sf potH, v4sf potL, double pot[]){
	pot[0] = (double)extract_ps(potH, 0)
		   + (double)extract_ps(potL, 0);
	pot[1] = (double)extract_ps(potH, 1)
		   + (double)extract_ps(potL, 1);
	pot[2] = (double)extract_ps(potH, 2)
		   + (double)extract_ps(potL, 2);
	pot[3] = (double)extract_ps(potH, 3)
		   + (double)extract_ps(potL, 3);
}

struct float2{
	float x, y;
};
static inline float2 float2_split(double x){
	float2 ret;
	x *= (1<<16);
	double xi = (int)x;
	double xf = x - xi;
	ret.x = xi * (1./(1<<16));
	ret.y = xf * (1./(1<<16));
	return ret;
}

struct Particle{
	float2 pos[3];
	float mass;
	float pad;

	Particle(double x[3], double m){
		pos[0] = float2_split(x[0]);
		pos[1] = float2_split(x[1]);
		pos[2] = float2_split(x[2]);
		mass = (float)m;
	}
	Particle(){
		pos[0].x = pos[0].y = pos[1].x = pos[1].y = pos[2].x = pos[2].y = mass = pad = 0.f;
	}
};

void gpupot(
        int rank,
        int istart,
        int ni,
		int n,
		double m[],
		double x[][3],
		double pot[]){
	double t0 = get_wtime();

	Particle *ptcl = new Particle[n+4];
	for(int i=0; i<n; i++){
		ptcl[i] = Particle(x[i], m[i]);
	}

    const int ibegin = istart - 1;
    const int iend = ni + istart - 1;

#pragma omp parallel for
	for(int i=ibegin; i<iend; i+=4){
		v4sf potH = _mm_setzero_ps();
		v4sf potL = _mm_setzero_ps();
		Particle *p = ptcl + i;
		v4sf xiH = _mm_set_ps(p[3].pos[0].x, p[2].pos[0].x, p[1].pos[0].x, p[0].pos[0].x);
		v4sf yiH = _mm_set_ps(p[3].pos[1].x, p[2].pos[1].x, p[1].pos[1].x, p[0].pos[1].x);
		v4sf ziH = _mm_set_ps(p[3].pos[2].x, p[2].pos[2].x, p[1].pos[2].x, p[0].pos[2].x);
		v4sf xiL = _mm_set_ps(p[3].pos[0].y, p[2].pos[0].y, p[1].pos[0].y, p[0].pos[0].y);
		v4sf yiL = _mm_set_ps(p[3].pos[1].y, p[2].pos[1].y, p[1].pos[1].y, p[0].pos[1].y);
		v4sf ziL = _mm_set_ps(p[3].pos[2].y, p[2].pos[2].y, p[1].pos[2].y, p[0].pos[2].y);
		for(int j=0; j<n; j++){
			v4sf jp0 = ((v4sf *)&ptcl[j])[0];
			v4sf jp1 = ((v4sf *)&ptcl[j])[1];
			v4sf xjH = _mm_shuffle_ps(jp0, jp0, 0x00);
			v4sf xjL = _mm_shuffle_ps(jp0, jp0, 0x55);
			v4sf yjH = _mm_shuffle_ps(jp0, jp0, 0xaa);
			v4sf yjL = _mm_shuffle_ps(jp0, jp0, 0xff);
			v4sf zjH = _mm_shuffle_ps(jp1, jp1, 0x00);
			v4sf zjL = _mm_shuffle_ps(jp1, jp1, 0x55);
			v4sf mj  = _mm_shuffle_ps(jp1, jp1, 0xaa);
			
			v4sf dx = (xjH - xiH) + (xjL - xiL);
			v4sf dy = (yjH - yiH) + (yjL - yiL);
			v4sf dz = (zjH - ziH) + (zjL - ziL);
			v4sf r2 = dx*dx + dy*dy + dz*dz;
			v4sf mask = _mm_cmplt_ps(_mm_setzero_ps(), r2);
			// Add small epsilon to avoid division by zero in rsqrt
			v4sf r2_safe = _mm_max_ps(r2, _mm_set1_ps(1.0e-30f));
			v4sf rinv = v4sf_rsqrt(r2_safe);
			rinv = _mm_and_ps(rinv, mask);
			rinv = rinv * mj;

			v4sf tmp = potH;
			potH = potH + rinv;
			potL = potL - ((potH - tmp) - rinv);
		}
		pot_reduce(potH, potL, pot+i);
	}

	delete [] ptcl;

	double t1 = get_wtime();
#ifdef PROFILE
	fprintf(stderr, "[R.%d AVX Pot.A] Ni %d  NTOT %d  pot(s) %f\n", rank,ni,n,t1 - t0);
#endif
    
}

extern "C"{
	void gpupot_(
            int *irank,
            int *istart,
            int *ni,
			int *n,
			double m[],
			double x[][3],
			double pot[]){
      gpupot(*irank,*istart,*ni,*n, m, x, pot);
	}
}
