#pragma once
#include <cassert>
#include <sys/time.h>
#define NAN_CHECK(value) assert((value) == (value))
#define PROFILE
static double get_wtime() { timeval tv; gettimeofday(&tv, nullptr); return tv.tv_sec + 1.e-6 * tv.tv_usec; }

// The float2 ("double-single") helpers below are used by gpupot.gpu.hip.  They
// exist in include/cuda_share.h but were dropped when this HIP header was
// generated, and HIPFLAGS does not add -I../include, so they are restored here.
// float2_split is called from the host constructor of Particle as well as from
// device code, hence __host__ __device__.
static __host__ __device__ float2 float2_split(double x){
	const int shift = 20;
	float2 ret;
	x *= (1<<shift);
	double xi = (int)x;
	double xf = x - xi;
	ret.x = xi * (1./(1<<shift));
	ret.y = xf * (1./(1<<shift));
	return ret;
}

static __device__ float2 float2_accum(float2 acc, float x){
	float tmp = acc.x + x;
	acc.y -= (tmp - acc.x) - x;
	acc.x = tmp;
	return acc;
}

static __device__ float2 float2_regularize(float2 acc){
	float tmp = acc.x + acc.y;
	acc.y = acc.y - (tmp - acc.x);
	acc.x = tmp;
	return acc;
}

static __device__ float2 float2_add(float2 a, float2 b){
	float tmp = a.x + b.x;
	a.y -= (tmp - a.x) - b.x - b.y;
	a.x = tmp;
	return a;
}
