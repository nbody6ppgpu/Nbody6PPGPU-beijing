#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([Nbody6++GPU-beijing], [Nov2024-hip], [spurzem@nao.cas.cn, kaiwu.astro@gmail.com])

m4_include([macro/ax_check_compile_flag.m4])
m4_include([macro/ax_fc_check_bounds.m4])

#AC_CONFIG_SRCDIR([helper_cuda.h.hip])

RESULT=
EXTRASRC=
EXTRAOBJ=

#Features controls
AC_ARG_ENABLE([mpi],
     [AS_HELP_STRING([--disable-mpi],
         [disable MPI support])],
     [],
     [enable_mpi=yes])

AC_ARG_ENABLE([mcmodel],
     [AS_HELP_STRING([--enable-mcmodel=ARG],
         [Use compiler option -mcmodel, avaiable arguments --> small; -->large; -->large
          (defalut is large)])],
     [ case "${enableval}" in
       large) enable_mcmodel=large ;;
       medium) enable_mcmodel=medium ;;
       small) enable_mcmodel=small ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-mcmodel) ;;
       esac],
     [enable_mcmodel=large])

AC_ARG_ENABLE([debug],
     [AS_HELP_STRING([--enable-debug],
         [Switch on debugging mode])],
     [],
     [enable_debug=no])

AC_ARG_ENABLE([gpu],
     [AS_HELP_STRING([--disable-gpu],
         [Disable GPU (HIP) acceleration support for regular force and potential])],
     [],
     [enable_gpu=yes])

AC_ARG_ENABLE([simd],
     [AS_HELP_STRING([--enable-simd=ARG],
         [Intel CPU SIMD instruction set acceleration for force and potential, avaiable arguments --> sse: Streaming SIMD Extensions; --> avx: Advanced Vector Extensions; --> no: switch off SIMD; (default is autodetecting)])],
     [ case "${enableval}" in
       avx) enable_simd=avx ;;
       sse) enable_simd=sse ;;
       no) enable_simd=no  ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-simd) ;;
       esac],
     [enable_simd=test])

AC_ARG_ENABLE([hdf5],
     [AS_HELP_STRING([--enable-hdf5],
         [Enable HDF5 output format support. (Ref: Cai, Maxwell, et al., 2015, submitted)])],
     [],
     [enable_hdf5=no])

AC_ARG_ENABLE([tt],
     [AS_HELP_STRING([--enable-tt],
         [Enable tidal tensor and analytical potential extension imported from NBODY6TT. (Ref: Renaud, F., et al., 2011, MNRAS, 418, 759; Renaud, F. & Gieles, M., 2015, MNRAS, 448, 3416)])],
     [],
     [enable_tt=no])

# Checks for programs.
# Fortran part

AC_LANG([Fortran])

local_fc=$FC
local_fcflags=$FCFLAGS
local_fflags=$FFLAGS

AC_PROG_FC
fc_check=$FC
FC=$local_fc

# Checks for MPI and HDF5 support
AS_IF([test "x$enable_mpi" != xno],
     [AS_IF([test "x$enable_hdf5" != xno],
       [AC_CHECK_PROG([FC], [h5pfc], [h5pfc], [none])],
       [AC_CHECK_PROG([FC], [mpif77], [mpif77], [none])])],
     [AS_IF([test "x$enable_hdf5" != xno],
      [AC_CHECK_PROG([FC], [h5fc], [h5fc], [none])],
      FC=$fc_check)
      local_fcflags=$local_fcflags' -I../extra_inc/nompi'])

AS_IF([test "x$enable_debug" != xno],
     [FCFLAGS=$local_fcflags' -g'
      AX_FC_CHECK_BOUNDS
      AX_CHECK_COMPILE_FLAG([-fbacktrace],[FCFLAGS=$FCFLAGS" -fbacktrace"])
      AX_CHECK_COMPILE_FLAG([-Wall],[FCFLAGS=$FCFLAGS" -Wall"])],
     [FCFLAGS=$local_fcflags" -O3"])

AX_CHECK_COMPILE_FLAG([-fPIC],
     [FCFLAGS=$FCFLAGS" -fPIC"],
     [AX_CHECK_COMPILE_FLAG([-fpic],
         [FCFLAGS=$FCFLAGS" -fpic"],
         [AC_MSG_FAILURE(
          [compiler option -fPIC/-fpic not avaiable, please use GNU or Intel Fortran compiler])]
         )])

AX_CHECK_COMPILE_FLAG([-mcmodel=$enable_mcmodel],
   [FCFLAGS=$FCFLAGS" -mcmodel="$enable_mcmodel],
   [AC_MSG_FAILURE(
       [$FC compiler option -mcmodel=$enable_mcmodel is not avaiable])])


AC_OPENMP
FCFLAGS=$FCFLAGS' '$OPENMP_FCFLAGS

FCFLAGS=$FCFLAGS' -I../include'

# Check FC
AS_IF([test "x$FC" == xnone],
     [AS_IF([test "x$enable_mpi" != xno],
        [AS_IF([test "x$enable_hdf5" != xno],
          [AC_MSG_FAILURE(
             [h5pfc not found, please install it and check environment PATH])],
          [AC_MSG_FAILURE(
            [mpif77 not found, please intall it and check environment PATH])])],
        [AS_IF([test "x$enable_hdf5" != xno],
          [AC_MSG_FAILURE(
             [h5fc not found, please install it and check environment PATH])],
          [AC_MSG_FAILURE(
            [Fortran compiler not found, please intall it and check environment PATH])])])])


#C++ part
AC_LANG([C++])
local_cxxflags=$CXXFLAGS
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_LIB([stdc++],[main],[LIBS=$LIBS' -lstdc++'])

CXXFLAGS=$local_cxxflags" -O3"

AC_OPENMP
CXXFLAGS=$CXXFLAGS' '$OPENMP_CXXFLAGS

CXXFLAGS=$CXXFLAGS' -I../include'

AX_CHECK_COMPILE_FLAG([-march=native],
    [CXXFLAGS=$CXXFLAGS" -march=native"],[])

	 
AX_CHECK_COMPILE_FLAG([-fPIC],
     CXXFLAGS=$CXXFLAGS" -fPIC",
     [AX_CHECK_COMPILE_FLAG([-fpic],
         CXXFLAGS=$CXXFLAGS" -fpic",
         [AC_MSG_FAILURE(
          [compiler option -fPIC/-fpic not avaiable, please use GNU or Intel C++ compiler])]
         )])

AX_CHECK_COMPILE_FLAG([-mcmodel=$enable_mcmodel],
   [CXXFLAGS=$CXXFLAGS" -mcmodel="$enable_mcmodel],
   [AC_MSG_FAILURE(
       [$CXX compiler option -mcmodel=$enable_mcmodel is not avaiable])])

# SIMD part
AS_IF([test "x$OPENMP_CXXFLAGS" != x],
   [AS_IF([test "x$enable_simd" == xtest],
       [c_simd_check="avx"]
       [AC_CHECK_FILE(["/proc/cpuinfo"],
         [AC_PROG_GREP
          SIMD_CHECK=`grep -c $c_simd_check /proc/cpuinfo`],
          SIMD_CHECK=-1)
         AC_MSG_CHECKING([whether current CPU support $c_simd_check])
         AS_IF([test "$SIMD_CHECK" -lt 0],
             [AC_MSG_RESULT([unknown])],
             [test "$SIMD_CHECK" -gt 0],
             [AC_MSG_RESULT([yes])],
             [AC_MSG_RESULT([no])])
         AX_CHECK_COMPILE_FLAG([-m$c_simd_check],
           [CXXFLAGS=$CXXFLAGS" -m"$c_simd_check
            RESULT=$RESULT"."$c_simd_check
	    enable_simd=$c_simd_check
           ],
           [c_simd_check="sse"]
           [AC_CHECK_FILE(["/proc/cpuinfo"],
              [AC_PROG_GREP
               SIMD_CHECK=`grep -c $c_simd_check /proc/cpuinfo`],
               SIMD_CHECK=-1)
            AC_MSG_CHECKING([whether current CPU support $c_simd_check])
            AS_IF([test "$SIMD_CHECK" -lt 0],
                  [AC_MSG_RESULT([unknown])],
                  [test "$SIMD_CHECK" -gt 0],
                  [AC_MSG_RESULT([yes])],
                  [AC_MSG_RESULT([no])])
            AX_CHECK_COMPILE_FLAG([-m$c_simd_check],
              [CXXFLAGS=$CXXFLAGS" -m"$c_simd_check
               RESULT=$RESULT"."$c_simd_check
	       enable_simd=$c_simd_check
              ],
              [enable_simd=no])])],
        [AS_IF([test "x$enable_simd" != xno],
               [AC_CHECK_FILE(["/proc/cpuinfo"],
                  [AC_PROG_GREP
                   SIMD_CHECK=`grep -c $enable_simd /proc/cpuinfo`],
                   SIMD_CHECK=-1)
                AC_MSG_CHECKING([whether current CPU support $enable_simd])
                AS_IF([test "$SIMD_CHECK" -lt 0],
                      [AC_MSG_RESULT([unknown])],
                      [test "$SIMD_CHECK" -gt 0],
                      [AC_MSG_RESULT([yes])],
                      [AC_MSG_RESULT([no])])
                AX_CHECK_COMPILE_FLAG([-m$enable_simd],
                  [CXXFLAGS=$CXXFLAGS" -m"$enable_simd
                   RESULT=$RESULT"."$enable_simd
                  ],
                  [AC_MSG_FAILURE(
                     [$CXX compiler option -m$enable_simd is not avaiable])])])])],
     enable_simd=no
     SIMD_CHECK=-2)

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday memset sqrt strchr])

#HIP part
# Checks for MPI support
AC_ARG_WITH([hip], 
     [AS_HELP_STRING([--with-hip=PREFIX],
         [Prefix of your HIP installation])],
     [hip_prefix=$withval],
     [hip_prefix="/opt/rocm/hip"])

AC_ARG_WITH([hip_sdk],
     [AS_HELP_STRING([--with-hip-sdk=PREFIX],
         [Prefix of your HIP SDKinstallation])],
     [hip_sdk_prefix=$withval],
     [hip_sdk_prefix="/opt/rocm/hip"])

AS_IF([test "x$enable_gpu" != xno],
     [AC_CHECK_PROG([HIPCC], [hipcc], [hipcc], [none])
      AS_IF([test "x$HIPCC" == xnone], 
            [AC_CHECK_PROG([HIPCC], [hipcc], [$hip_prefix/bin/hipcc], [none], [$hip_prefix/bin])
             AS_IF([test "x$HIPCC" == xnone], 
                   [AC_MSG_FAILURE(
                       [Cannot find HIP compiler hipcc, please check whether hipcc is in environment PATH or use --with-hip=PREFIX to provide the PATH of HIP installation])])],
            [AC_PROG_SED
             hip_prefix=`hipconfig --rocmpath`])
      AC_CHECK_LIB([amdhip64], [main], 
          [LIBS=$LIBS' -lamdhip64'],
          [AC_CHECK_FILE(["$hip_prefix/lib/libamdhip64.so"],
              LIBS=$LIBS" -L$hip_prefix/lib -lamdhip64 -L$hip_prefix/llvm/lib",
              [AC_MSG_FAILURE(
                  [Cannot find HIP library -lamdhip64, please provide correct hip PREFIX by using --with-hip=PREFIX])])])
      HIPFLAGS=" -I ../extra_inc/cuda"
      RESULT=$RESULT".hip"
      AC_PROG_AWK
    #   AC_MSG_CHECKING([for HIP version])
    #   HIP_VERSION=`$HIPCC --version | sed -n 's/^HIP version: //p'`
    #   HIP_VERSION_I=`echo $HIP_VERSION|awk -F '.' '{print $1}'`

      FCFLAGS=$FCFLAGS' $(GPU_FLAGS)'
      EXTRASRC=$EXTRASRC' $(GPU_FSOURCES)'
      AS_IF([test "x$OPENMP_CXXFLAGS" != x ],
           EXTRAOBJ=$EXTRAOBJ' $(HIP_OBJECTS)', 
           EXTRAOBJ=$EXTRAOBJ' $(HIP_OBJECTS_SINGLE)')
      LDFLAGS=$LDFLAGS' -L'$HIP_PREFIX
      dnl AX_CHECK_COMPILE_FLAG([-Xcompiler],[],
      dnl    [AC_MSG_FAILURE(
      dnl       [Current hipcc do not support -Xcompiler])])
     ])

#FLAGS
# MPI FLAGS
AS_IF([test "x$enable_mpi" != xno],
      FCFLAGS=$FCFLAGS' $(MPI_FLAGS)'
      # GCC >= 10 requires -fallow-argument-mismatch flag in MPI
      # https://gcc.gnu.org/gcc-10/changes.html#fortran
      AX_CHECK_COMPILE_FLAG([-fallow-argument-mismatch], [MPI_FLAGS=$MPI_FLAGS" -fallow-argument-mismatch"])
      EXTRASRC=$EXTRASRC' $(MPI_FSOURCES)'
      RESULT=$RESULT'.mpi')

#SIMD
AS_IF([test "x$enable_simd" != xno],
      [FCFLAGS=$FCFLAGS' ${SIMD_FLAGS}']
      [CXXCOMP=`echo x$CXX|grep -c icc`]
      [AS_IF([test $CXXCOMP -gt 0],
         CXXFLAGS=$CXXFLAGS' $(INTEL_FLAGS)')])

#GPU
AS_IF([test "x$enable_simd" == xavx],
      [AS_IF([test "x$enable_gpu" != xno],
             EXTRAOBJ=$EXTRAOBJ' ${IRRAVX_OBJECTS}',
             EXTRAOBJ=$EXTRAOBJ' ${AVX_OBJECTS}'
             EXTRASRC=$EXTRASRC' $(GPU_FSOURCES)'
             FCFLAGS=$FCFLAGS' $(GPU_FLAGS)')],
      [test "x$enable_simd" == xsse],
      [AS_IF([test "x$enable_gpu" != xno],
             EXTRAOBJ=$EXTRAOBJ' ${IRRSSE_OBJECTS}',
             EXTRAOBJ=$EXTRAOBJ' ${SSE_OBJECTS}'
             EXTRASRC=$EXTRASRC' $(GPU_FSOURCES)'
             FCFLAGS=$FCFLAGS' $(GPU_FLAGS)')])

#OMP
AS_IF([test "x$OPENMP_FCFLAGS" != x ],
      FCFLAGS=$FCFLAGS' ${OMP_FLAGS}',
      FCFLAGS=$FCFLAGS' -I../extra_inc/noomp')

AS_IF([test "x$OPENMP_CXXFLAGS" != x ],
      CXXFLAGS=$CXXFLAGS' ${OMP_FLAGS}',
      CXXFLAGS=$CXXFLAGS' -I ../extra_inc/noomp')

# HDF5 format
AS_IF([test "x$enable_hdf5" != xno],
       [FCFLAGS=$FCFLAGS' $(HDF5_FLAGS)'
        EXTRASRC=$EXTRASRC' $(HDF5_FSOURCES)'
        RESULT=$RESULT'.hdf5'])

# TT support
AS_IF([test "x$enable_tt" != xno],
       [FCFLAGS=$FCFLAGS' $(TT_FLAGS)'
        EXTRASRC=$EXTRASRC' $(TT_FSOURCES)'
	RESULT=$RESULT'.tt'])

# DEBUG FLAGS
AS_IF([test "x$enable_debug" != xno],
      FCFLAGS=$FCFLAGS' $(DEBUG_FLAGS)'
      RESULT=$RESULT'.db')


# params.h
AC_ARG_WITH([par], 
     [AS_HELP_STRING([--with-par=size],
         [Size of basic parameters (NMAX, KMAX, LMAX, MMAX) in params.h. Avaiable values: 10m, 8m, 4m, 1m, 512k, 256k, 128k, 64k, 32k, 16k, 8k, custom. The value of this size is NMAX, with prefix 'b' the K.S. binary size KMAX is half of NMAX and NMAX is set to 1.5*NMAX. LMAX is 600 above 64k, 500 above 16k. MMAX is 1024 above 64k, 512 above 16k. MMAX is doubled with prefix 'b'. Defaultsize is b1m])],
     [par_size=$withval],
     [par_size=b1m])

AS_IF([test "x$par_size" == x16m],
      [NMAX=16777216
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb16m],
      [NMAX=25165824
       KMAX=8388608
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x10m],
      [NMAX=11048576
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb10m],
      [NMAX=16572864
       KMAX=5524288
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x8m],
      [NMAX=8388608
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb8m],
      [NMAX=12582912
       KMAX=4194304
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x4m],
      [NMAX=4194304
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb4m],
      [NMAX=6291456
       KMAX=2097152
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x1m],
      [NMAX=1048576
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb1m],
      [NMAX=1572864
       KMAX=524288
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x512k],
      [NMAX=524288
       KMAX=65536
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb512k],
      [NMAX=786432
       KMAX=262144
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x256k],
      [NMAX=262144
       KMAX=65536
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb256k],
      [NMAX=393216
       KMAX=131072
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x128k],
      [NMAX=131072
       KMAX=65536
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb128k],
      [NMAX=196608
       KMAX=65536
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x64k],
      [NMAX=65536
       KMAX=8192
       LMAX=500
       MMAX=512],
      [test "x$par_size" == xb64k],
      [NMAX=98304
       KMAX=32768
       LMAX=500
       MMAX=1024],
      [test "x$par_size" == x32k],
      [NMAX=32768
       KMAX=8192
       LMAX=500
       MMAX=512],
      [test "x$par_size" == xb32k],
      [NMAX=49152
       KMAX=16384
       LMAX=500
       MMAX=1024],
      [test "x$par_size" == x16k],
      [NMAX=16384
       KMAX=2048
       LMAX=400
       MMAX=256],
      [test "x$par_size" == xb16k],
      [NMAX=24576
       KMAX=8192
       LMAX=400
       MMAX=512],
      [test "x$par_size" == x8k],
      [NMAX=8192
       KMAX=2048
       LMAX=400
       MMAX=256],
      [test "x$par_size" == xb8k],
      [NMAX=12288
       KMAX=4096
       LMAX=400
       MMAX=512],
      [test "x$par_size" == xcustom],
      [flag_input_par=yes],
      [NMAX=24576
       KMAX=8192
       LMAX=400
       MMAX=512])

AC_ARG_WITH([nmax], 
     [AS_HELP_STRING([--with-nmax=size],
         [Set Maximum number of particles (will overwrite value set in --with-par)])],
     [NMAX=$withval],
     [flag_par_miss=yes])

AC_ARG_WITH([kmax], 
     [AS_HELP_STRING([--with-kmax=size],
         [Set Maximum number of K.S. particles (will overwrite value set in --with-par)])],
     [KMAX=$withval],
     [flag_par_miss=yes])

AC_ARG_WITH([lmax], 
     [AS_HELP_STRING([--with-lmax=size],
         [Set Maximum number of neighbor list particles (will overwrite value set in --with-par)])],
     [LMAX=$withval],
     [flag_par_miss=yes])

AC_ARG_WITH([mmax], 
     [AS_HELP_STRING([--with-mmax=size],
         [Set Maximum number of merger (stable triple) particles (will overwrite value set in --with-par)])],
     [MMAX=$withval],
     [flag_par_miss=yes])

AS_IF([test "x$flag_input_par" == xyes],
  [AS_IF([test "x$flag_par_miss" == xyes],
        [AC_MSG_FAILURE(
             [Parameters are not completely set, please make sure that all four option (--with-nmax, --with-kmax, --with-lmax, --with-mmax) are set])])])

AC_SUBST([HIPFLAGS])	
AC_SUBST([FCFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([RESULT])
AC_SUBST([EXTRASRC])
AC_SUBST([EXTRAOBJ])
AC_SUBST([EXTRARESULT])
AC_SUBST([MPI_FLAGS])

AC_SUBST([NMAX])
AC_SUBST([KMAX])
AC_SUBST([LMAX])
AC_SUBST([MMAX])

AC_MSG_NOTICE([----------------------Configure Summary--------------------])
AC_MSG_NOTICE([--Features:])
AC_MSG_NOTICE([     Using MPI:      $enable_mpi])
AC_MSG_NOTICE([     Using GPU:      $enable_gpu])
AC_MSG_NOTICE([     HIP version:   $HIP_VERSION])
AC_MSG_NOTICE([     Using SIMD:     $enable_simd])
AS_IF([test "x$enable_simd" != xno],
  [AS_IF([test "$SIMD_CHECK" == -1],
      [AC_MSG_WARN([      No /proc/cpuinfo found, cannot comfirm whether current CPU support SIMD])],
        [test "$SIMD_CHECK" == 0],
      [AC_MSG_WARN([      SIMD option used now is not supported in current CPU])])],
        [test "$SIMD_CHECK" == -2],
      [AC_MSG_WARN([      SIMD is suppressed since OpenMP is switched off])])
AC_MSG_NOTICE([     mcmodel:        $enable_mcmodel])
AC_MSG_NOTICE([     OpenMP option:  $OPENMP_FCFLAGS])
AC_MSG_NOTICE([     Debug mode:     $enable_debug])
AC_MSG_NOTICE([     HDF5 format:    $enable_hdf5])
AC_MSG_NOTICE([     TT package:     $enable_tt])
AC_MSG_NOTICE([--Compilers:])
AC_MSG_NOTICE([     C++ compiler:      $CXX])
AC_MSG_NOTICE([     Fortran compiler:  $FC])
AC_MSG_NOTICE([     HIP compiler:      $HIPCC])
AC_MSG_NOTICE([--Targets:])
AC_MSG_NOTICE([     Program name:   nbody6++$RESULT])
AC_MSG_NOTICE([     Install path:   $prefix])
AC_MSG_NOTICE([--Parameters:])
AC_MSG_NOTICE([     NMAX:           $NMAX])
AC_MSG_NOTICE([     KMAX:           $KMAX])
AC_MSG_NOTICE([     LMAX:           $LMAX])
AC_MSG_NOTICE([     MMAX:           $MMAX])
AC_MSG_NOTICE([-----------------------------------------------------------])
AC_MSG_NOTICE([Reference paper: Wang, L., et al., 2015, MNRAS, 450, 4070 (http://arxiv.org/abs/1504.03687)])
AC_MSG_NOTICE([Email: Rainer Spurzem<spurzem@nao.cas.cn>])
AC_MSG_NOTICE([Check updates or ask anything in https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing])
AC_MSG_NOTICE([-----------------------------------------------------------])

AC_CONFIG_FILES([Makefile build/Makefile include/params.h])
AC_OUTPUT
