#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include "simd_define.h"

/* Include SSE intrinsics if using SIMDe or native x86 */
#ifdef NBODY_USE_SIMDE
/* SIMDe headers already included via simd_define.h */
#elif !defined(__USE_INTEL)
#include <xmmintrin.h>  /* SSE */
#include <emmintrin.h>  /* SSE2 */
#include <pmmintrin.h>  /* SSE3 */
#endif

#define TMAX 1024 // maximum number of threads
#if 1
#include <omp.h>
#else
static inline int omp_get_num_threads(){return 1;}
static inline int omp_get_thread_num() {return 0;}
#endif

#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
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

static double time_send, time_grav;
static long long numInter;
static int isend,icall,ini;

template <class T>
struct myvector{
	int num;
	T *val;
	myvector(){
		num = 0;
		val = NULL;
	}
	~myvector(){
		delete [] val;
	}
	void clear(){
		num = 0;
	}
	void reserve(size_t count){
		val = new T[count];
	}
	void free(){
		delete [] val;
	}
	void push_back(const T &t){
		val[num++] = t;
	}
	size_t size(){
		return num;
	}
	T &operator[](int i){
		return val[i];
	}
};

//typedef float v4sf __attribute__ ((vector_size(16)));
static inline v4sf v4sf_rsqrt(v4sf x){
	v4sf y = _mm_rsqrt_ps(x);
	return (_mm_set1_ps(-0.5f) * y) * 
			(x*y*y + _mm_set1_ps(-3.f));
}
// #include "v4sf.h"

struct Jparticle{
	float x[3];
	float m;
	float v[3];
	float pad;
	Jparticle() {}
	Jparticle(double mj, double xj[3], double vj[3]){
		x[0] = xj[0];
		x[1] = xj[1];
		x[2] = xj[2];
		m    = mj;
		v[0] = vj[0];
		v[1] = vj[1];
		v[2] = vj[2];
	}
};
static Jparticle *jp_host;
// static int *nblist;
// static int *nblistbuf[TMAX];
static myvector<int> nblist[TMAX][4];
static int nbody, nbodymax;

void GPUNB_open(int nbmax, int irank){
	// std::cout << "Open GPUNB " << nbmax << std::endl;
	time_send = time_grav = 0.0;
	numInter = 0;
	isend = icall = ini = 0;
	nbodymax = nbmax;
	jp_host = new Jparticle[nbmax+3];
    int numCPU = -1;
#pragma omp parallel
	{
		int nth = omp_get_num_threads();
		assert(nth <= TMAX);
		int tid = omp_get_thread_num();
		for(int v=0; v<4; v++){
			nblist[tid][v].reserve(nbmax);
		}
        if(tid == 0) numCPU = omp_get_num_threads();
	}
#ifdef PROFILE
    fprintf(stderr, "# Open SSE regular force - rank: %d; threads: %d\n", irank, numCPU);
#endif
}

void GPUNB_close(){
	// std::cout << "Close GPUNB" << std::endl;
	delete [] jp_host;
#pragma omp parallel
	{
		// int tid = omp_get_thread_num();
		// for(int v=0; v<4; v++){
		// 	nblist[tid][v]::~vector();
		// }
		/*
		int tid = omp_get_thread_num();
		for(int v=0; v<4; v++){
			nblist[tid][v].free();
		}*/
	}
	nbodymax = 0;

// #ifdef PROFILE
// 	std::cerr << "***********************" << std::endl;
// 	std::cerr << "time send : " << time_send << " sec " << std::endl;
// 	std::cerr << "time grav : " << time_grav << " sec " << std::endl;
// 	std::cerr << 60.e-9 * numInter / time_grav << " Gflops (gravity part only)" << std::endl;
// 	std::cerr << "***********************" << std::endl;
// #endif
}

void GPUNB_send(
		int nj,
		double mj[],
		double xj[][3],
		double vj[][3]){
	time_send -= get_wtime();
	nbody = nj;
    ::isend++;
	// std::cout << "gpu send: " << nbody << " " << nbodymax << std::endl;
	assert(nbody <= nbodymax);
#pragma omp parallel for
	for(int j=0; j<nj; j++){
		jp_host[j] = Jparticle(mj[j], xj[j], vj[j]);
	}
	time_send += get_wtime();
}

void GPUNB_regf(
		int ni,
		double h2d[],
		double dtr[],
		double xid[][3],
		double vid[][3],
		double acc[][3],
		double jrk[][3],
		double pot[],
		int lmax,
		int nbmax,
		int *listbase,
        int m_flag){
	// std::cout << " Call GPUNB_regf " << ni << std::endl;
	time_grav -= get_wtime();
	numInter += ni * nbody;
    ::icall++;
    ::ini +=ni;
#pragma omp parallel for
	for(int i=0; i<ni; i+=4){
		int tid = omp_get_thread_num();
		nblist[tid][0].clear();
		nblist[tid][1].clear();
		nblist[tid][2].clear();
		nblist[tid][3].clear();
		int nii = std::min(4, ni-i);

		v4sf xi  = _mm_set_ps((float)xid[i+3][0], (float)xid[i+2][0], (float)xid[i+1][0], (float)xid[i+0][0]); 
		v4sf yi  = _mm_set_ps((float)xid[i+3][1], (float)xid[i+2][1], (float)xid[i+1][1], (float)xid[i+0][1]); 
		v4sf zi  = _mm_set_ps((float)xid[i+3][2], (float)xid[i+2][2], (float)xid[i+1][2], (float)xid[i+0][2]); 
		v4sf vxi = _mm_set_ps((float)vid[i+3][0], (float)vid[i+2][0], (float)vid[i+1][0], (float)vid[i+0][0]); 
		v4sf vyi = _mm_set_ps((float)vid[i+3][1], (float)vid[i+2][1], (float)vid[i+1][1], (float)vid[i+0][1]); 
		v4sf vzi = _mm_set_ps((float)vid[i+3][2], (float)vid[i+2][2], (float)vid[i+1][2], (float)vid[i+0][2]); 
		v4sf h2i = _mm_set_ps((float)h2d[i+3], (float)h2d[i+2], (float)h2d[i+1], (float)h2d[i+0]); 
		/* Mask arrays for partial vector operations */
		v4sf h2mask_val;
		switch(nii) {
			case 0: h2mask_val = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f); break;
			case 1: h2mask_val = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f); break;
			case 2: h2mask_val = _mm_set_ps(0.0f, 0.0f, 1.0f, 1.0f); break;
			case 3: h2mask_val = _mm_set_ps(0.0f, 1.0f, 1.0f, 1.0f); break;
			default: h2mask_val = _mm_set_ps(1.0f, 1.0f, 1.0f, 1.0f); break;
		}
		h2i = h2i * h2mask_val;
		v4sf dtri = _mm_set_ps((float)dtr[i+3], (float)dtr[i+2], (float)dtr[i+1], (float)dtr[i+0]); 
		v4sf Ax = _mm_setzero_ps();
		v4sf Ay = _mm_setzero_ps();
		v4sf Az = _mm_setzero_ps();
		v4sf Jx = _mm_setzero_ps();
		v4sf Jy = _mm_setzero_ps();
		v4sf Jz = _mm_setzero_ps();
		v4sf poti = _mm_setzero_ps();
		v4sf *jpp = (v4sf *)jp_host;
		for(int j=0; j<nbody; j++, jpp+=2){
			v4sf jp0 = jpp[0];
			v4sf jp1 = jpp[1];

			v4sf xj = _mm_shuffle_ps(jp0, jp0, 0x00);
			v4sf yj = _mm_shuffle_ps(jp0, jp0, 0x55);
			v4sf zj = _mm_shuffle_ps(jp0, jp0, 0xaa);
			v4sf mj = _mm_shuffle_ps(jp0, jp0, 0xff);
			v4sf vxj = _mm_shuffle_ps(jp1, jp1, 0x00);
			v4sf vyj = _mm_shuffle_ps(jp1, jp1, 0x55);
			v4sf vzj = _mm_shuffle_ps(jp1, jp1, 0xaa);

			v4sf dx = xj - xi;
			v4sf dy = yj - yi;
			v4sf dz = zj - zi;
			v4sf dvx = vxj - vxi;
			v4sf dvy = vyj - vyi;
			v4sf dvz = vzj - vzi;

			v4sf dxp = dx + dtri * dvx;
			v4sf dyp = dy + dtri * dvy;
			v4sf dzp = dz + dtri * dvz;

			v4sf r2 = dx*dx + dy*dy + dz*dz;
			v4sf rv = dx*dvx + dy*dvy + dz*dvz;
			v4sf r2p = dxp*dxp + dyp*dyp + dzp*dzp;
            v4sf mask;
            //          v4sf mask = _mm_cmplt_ps(r2, h2i);
            if(m_flag) {
              v4sf mh2i = mj * h2i;
              mask = _mm_cmplt_ps(
                       _mm_min_ps(r2,r2p), mh2i);
            }
            else {
              mask = _mm_cmplt_ps(
                       _mm_min_ps(r2,r2p), h2i);
            }
			int bits = _mm_movemask_ps(mask);
			// mj = _mm_andnot_ps(mask, mj);
			if(bits){
				if (bits&1) nblist[tid][0].push_back(j);
				if (bits&2) nblist[tid][1].push_back(j);
				if (bits&4) nblist[tid][2].push_back(j);
				if (bits&8) nblist[tid][3].push_back(j);
			}

			// Add small epsilon to avoid division by zero in rsqrt
			v4sf r2_safe = _mm_max_ps(r2, _mm_set1_ps(1.0e-30f));
			v4sf rinv1 = v4sf_rsqrt(r2_safe);
			rinv1 = _mm_andnot_ps(mask, rinv1);
			// v4sf rinv1 = _mm_rsqrt_ps(r2);
			v4sf rinv2 = rinv1 * rinv1;
			rinv1 = rinv1 * mj;
			poti = poti + rinv1;
			v4sf rinv3 = rinv1 * rinv2;
			rv = rv * (_mm_set1_ps(-3.f) * rinv2);

			Ax = Ax + rinv3 * dx;
			Ay = Ay + rinv3 * dy;
			Az = Az + rinv3 * dz;
			Jx = Jx + rinv3 * (dvx + rv * dx);
			Jy = Jy + rinv3 * (dvy + rv * dy);
			Jz = Jz + rinv3 * (dvz + rv * dz);
		} // for(j)
		union {
			struct{
				__m128 Ax, Ay, Az, Jx, Jy, Jz, Pot;  /* Use underlying type for union compatibility */
			};
			struct{
				float acc[3][4], jrk[3][4], pot[4];
			};
		} u;
		u.Ax = Ax;
		u.Ay = Ay;
		u.Az = Az;
		u.Jx = Jx;
		u.Jy = Jy;
		u.Jz = Jz;
		u.Pot = poti;
		for(int ii=0; ii<nii; ii++){
			for(int k=0; k<3; k++){
				acc[i+ii][k] = u.acc[k][ii];
				jrk[i+ii][k] = u.jrk[k][ii];
			}
			pot[i+ii] = u.pot[ii];
			int nnb = nblist[tid][ii].size();
			int *nnbp = listbase + lmax * (i+ii);
			int *nblistp = nnbp + 1;
			if(nnb > nbmax){
				*nnbp = -nnb;
			}else{
				*nnbp = nnb;
				for(int k=0; k<nnb; k++){
					nblistp[k] = nblist[tid][ii][k];
				}
			}
		}
	}
	// printf("gpu: %e %e %e %d\n", xid[0][0], acc[0][0], jrk[0][0], *listbase);
#if 0
	if(ni > 0){
		FILE *fp = fopen("Force.sse", "w");
		assert(fp);
		for(int i=0; i<ni; i++){
			int nnb =  listbase[i*lmax];
			fprintf(fp, "%d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %d\n",
					i, acc[i][0], acc[i][1], acc[i][2], 
					   jrk[i][0], jrk[i][1], jrk[i][2], nnb);
		}
		fprintf(fp, "\n");
		fclose(fp);
		exit(1);
	}
#endif
	time_grav += get_wtime();
}

void GPUNB_profile(int irank) {
#ifdef PROFILE
  if(icall) {
    // R: rank; D: GPU device number;
    // Nsend: number of call gpunb_send between two profile check points (adjust time interval);
    // Ngrav: number of call gpunb_ref between two profile check points;
    // <Ni>:  averaged i particle number per regular block;
    // send: j particle sending time;
    // grav: force calculation time;
    // Perf: performance for gpu regular calculation
    int avg_ni = (isend > 0) ? ini/isend : 0;
    fprintf(stderr,"[R.%d SSE Reg.F ] Nsend %d  Ngrav %d  <Ni> %d   send(s) %f grav(s) %f  Perf.(Gflops) %f\n",irank,isend,icall,avg_ni,time_send,time_grav,60.e-9*numInter/time_grav);
  }
  time_send = time_grav = 0.0;
  numInter = 0;
  icall = ini = isend= 0;
#else
  return;
#endif
}

extern "C" {
    void gpunb_open_( int *nbmax, int *irank){
        GPUNB_open(*nbmax,*irank);
	}
	void gpunb_close_(){
		GPUNB_close();
	}
	void gpunb_send_(
			int *nj,
			double mj[],
			double xj[][3],
			double vj[][3]){
		GPUNB_send(*nj, mj, xj, vj);
	}
	void gpunb_regf_(
			int *ni,
			double h2[],
			double dtr[],
			double xi[][3],
			double vi[][3],
			double acc[][3],
			double jrk[][3],
			double pot[],
			int *lmax,
			int *nbmax,
			int *list,
            int *m_flag){ // list[][lmax]
        GPUNB_regf(*ni, h2, dtr, xi, vi, acc, jrk, pot, *lmax, *nbmax, list, *m_flag);
	}
    void gpunb_profile_(int *irank){
      GPUNB_profile(*irank);
    }
}
