#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([Nbody6++], [1.0], [longwang.astro@gmail.com])
#AC_CONFIG_MACRO_DIR([./macro])

# Detect host architecture for platform-specific settings
AC_CANONICAL_HOST

source_root=`cd "$srcdir" && pwd`
source_parent=`cd "$source_root/.." && pwd`
simde_default_path=${source_parent%/}/simde
arm_configure_example="./configure --with-par=b1m --with-simde=\"$simde_default_path\" --enable-simd=sse --disable-gpu"
x86_configure_example='./configure --with-par=b1m --enable-mcmodel=large --disable-gpu'

AS_CASE([$host_cpu],
  [x86_64|amd64], [target_architecture=x86_64],
  [aarch64|arm64], [target_architecture=ARM64],
  [AC_MSG_ERROR([unsupported target architecture '$host_cpu' (canonical host: $host).
Supported target architectures: x86_64 and ARM64/aarch64.
x86_64 example: $x86_configure_example
ARM64 example: $arm_configure_example])])

AC_MSG_NOTICE([target architecture: $target_architecture])
AC_MSG_NOTICE([canonical host: $host])
AS_IF([test "x$build" != "x$host"], [
  AC_MSG_NOTICE([cross compilation: build=$build, host=$host])
])

m4_include([macro/ax_check_compile_flag.m4])
m4_include([macro/ax_fc_check_bounds.m4])
#m4_include([macro/ax_check_cuda.m4])
#AC_CONFIG_SRCDIR([helper_cuda.h])

RESULT=
EXTRASRC=
EXTRAOBJ=

#Save environment
##AC_ARG_ENABLE([local-env],
##   [AS_HELP_STRING([--enable-local-env],
##        [Use local environment listed in last part instead of configure])],
##   [],
##   [enable_local_env=no])

##AS_IF([test "x$enable_local_env" != xno ],
##     [AS_IF([test "x$FC" != x], [local_fc=$FC])
##      AS_IF([test "x$FCFLAGS" != x], [local_fcflags=$FCFLAGS])
##      AS_IF([test "x$CXXFLAGS" != x], [local_cxxflags=$CXXFLAGS])],
##    [local_fc=
##     local_fcflags=
##     local_cxxflags=])

#Features controls
AC_ARG_ENABLE([mpi],
     [AS_HELP_STRING([--disable-mpi],
         [disable MPI support])],
     [],
     [enable_mpi=yes])

AC_ARG_ENABLE([mcmodel],
     [AS_HELP_STRING([--enable-mcmodel=ARG],
         [Use compiler option -mcmodel, avaiable arguments --> small; -->medium; -->large; -->no
          (default is large on x86, disabled on ARM64 where mcmodel is not supported)])],
     [ mcmodel_was_specified=yes
       case "${enableval}" in
       large) enable_mcmodel=large ;;
       medium) enable_mcmodel=medium ;;
       small) enable_mcmodel=small ;;
       no) enable_mcmodel=no ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-mcmodel) ;;
       esac],
     [mcmodel_was_specified=no
      AS_CASE([$target_architecture],
       [ARM64], [enable_mcmodel=no],
       [enable_mcmodel=large])])

AC_ARG_ENABLE([debug],
     [AS_HELP_STRING([--enable-debug],
         [Switch on debugging mode])],
     [],
     [enable_debug=no])

AC_ARG_ENABLE([gpu],
     [AS_HELP_STRING([--disable-gpu],
         [Disable CUDA/HIP acceleration support for regular force and potential])],
     [],
     [enable_gpu=yes])

AC_ARG_WITH([gpu-backend],
     [AS_HELP_STRING([--with-gpu-backend=auto|cuda|hip],
         [Select GPU backend (default: auto, preferring CUDA over HIP)])],
     [gpu_backend=$withval; gpu_backend_explicit=yes],
     [gpu_backend=auto; gpu_backend_explicit=no])
AS_CASE([$gpu_backend], [auto|cuda|hip], [],
   [AC_MSG_ERROR([bad value $gpu_backend for --with-gpu-backend])])
AS_IF([test "x$enable_gpu" = xno && test "x$gpu_backend_explicit" = xyes && test "x$gpu_backend" != xauto],
   [AC_MSG_ERROR([--disable-gpu conflicts with --with-gpu-backend=$gpu_backend])])

AC_ARG_ENABLE([simd],
     [AS_HELP_STRING([--enable-simd=ARG],
         [SIMD instruction set acceleration for force and potential, available arguments --> sse: Streaming SIMD Extensions (default on ARM64, via SIMDe); --> avx: Advanced Vector Extensions (x86_64 only); --> no: switch off SIMD; this also disables OpenMP, since the resulting plain nbint.F force routine is not thread-safe under OpenMP; (default is autodetecting on x86_64 and sse on ARM64)])],
     [ simd_was_specified=yes
       case "${enableval}" in
       avx) enable_simd=avx ;;
       sse) enable_simd=sse ;;
       no) enable_simd=no  ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-simd) ;;
       esac],
     [simd_was_specified=no
      AS_CASE([$target_architecture],
       [ARM64], [enable_simd=sse],
       [enable_simd=test])])

# SIMDe support for ARM and cross-platform SIMD
AC_ARG_WITH([simde],
     [AS_HELP_STRING([--with-simde@<:@=PATH@:>@],
         [Use SIMDe headers. On ARM64, the default PATH is the source sibling ../simde; a bare --with-simde uses system include paths; --without-simde is invalid unless --enable-simd=no is also given.])],
     [simde_was_specified=yes
      with_simde=$withval],
     [simde_was_specified=no
      AS_IF([test "x$target_architecture" = xARM64 && test "x$enable_simd" != xno],
        [with_simde=$simde_default_path],
        [with_simde=no])])

AS_IF([test "x$target_architecture" = xARM64], [
  AS_IF([test "x$enable_simd" = xavx],
    [AC_MSG_ERROR([option --enable-simd=avx conflicts with target architecture ARM64; AVX is x86_64-only.
ARM64 example: $arm_configure_example])])
  AS_IF([test "x$enable_mcmodel" != xno],
    [AC_MSG_ERROR([option --enable-mcmodel=$enable_mcmodel conflicts with target architecture ARM64; use --enable-mcmodel=no or omit the option.
ARM64 example: $arm_configure_example])])
  AS_IF([test "x$enable_simd" != xno && test "x$with_simde" = xno],
    [AC_MSG_ERROR([option --without-simde conflicts with target architecture ARM64 unless combined with --enable-simd=no; SIMDe is required for the SSE kernels.
ARM64 example: $arm_configure_example])])
])

# Normalize relative path for SIMDe to an absolute path (build runs in ./build)
AS_IF([test "x$with_simde" != xno && test "x$with_simde" != xyes && test "x$with_simde" != x],
   [AS_CASE([$with_simde],
      [/*], [],
      [*], [AS_IF([test -d "$with_simde"],
             [with_simde=`cd "$with_simde" && pwd`],
             [with_simde="$source_root/$with_simde"])])])

AC_ARG_ENABLE([hdf5],
     [AS_HELP_STRING([--disable-hdf5],
         [Disable HDF5 output format support. HDF5 is auto-detected and used by
default, since all current output formats are built on it and the legacy
ASCII output is unmaintained; disabling it is NOT recommended.
(Ref: Cai, Maxwell, et al., 2015, submitted)])],
     [],
     [enable_hdf5=yes])

AC_ARG_VAR([H5FC], [HDF5 Fortran wrapper used to discover the serial HDF5 installation])

AC_ARG_ENABLE([tt],
     [AS_HELP_STRING([--enable-tt],
         [Enable tidal tensor and analytical potential extension imported from NBODY6TT. (Ref: Renaud, F., et al., 2011, MNRAS, 418, 759; Renaud, F. & Gieles, M., 2015, MNRAS, 448, 3416)])],
     [],
     [enable_tt=no])

AC_ARG_ENABLE([tools],
    [AS_HELP_STRING([--enable-tools],
         [Compile extral tools for data analysis, see User_manual.pdf for detail])],
    [enable_tools=yes],
    [enable_tools=no])

# Checks for programs.
# Fortran part

AC_LANG([Fortran])

local_fc=$FC
local_fcflags=$FCFLAGS
local_fflags=$FFLAGS

AC_PROG_FC
fc_check=$FC
FC=$local_fc

# Check the Fortran compiler used by the simulation.  HDF5 is discovered
# separately below so an MPI build can link the serial HDF5 library.
AS_IF([test "x$enable_mpi" != xno],
     [AC_CHECK_PROG([FC], [mpif77], [mpif77], [none])],
     [FC=$fc_check
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

AS_IF([test "x$enable_mcmodel" != xno],
   [AX_CHECK_COMPILE_FLAG([-mcmodel=$enable_mcmodel],
      [FCFLAGS=$FCFLAGS" -mcmodel="$enable_mcmodel],
      [AC_MSG_FAILURE(
          [$FC compiler option -mcmodel=$enable_mcmodel is not avaiable])])])


AC_OPENMP

FCFLAGS=$FCFLAGS' -I../include'

# Check FC
AS_IF([test "x$FC" == xnone],
     [AS_IF([test "x$enable_mpi" != xno],
        [AC_MSG_FAILURE(
           [mpif77 not found, please install it and check environment PATH])],
        [AC_MSG_FAILURE(
           [Fortran compiler not found, please install it and check environment PATH])])])

HDF5_FCFLAGS=
HDF5_LIBS=
HDF5_PARALLEL=disabled

AS_IF([test "x$enable_hdf5" != xno], [
  AC_PATH_PROG([H5FC], [h5fc], [none])
  AS_IF([test "x$H5FC" = xnone],
    [AC_MSG_FAILURE(
       [h5fc not found. HDF5 is required by default: install the serial HDF5
Fortran development package (e.g. libhdf5-dev on Debian/Ubuntu), load a
serial HDF5 module, or set H5FC=/path/to/h5fc. Alternatively, reconfigure
with --disable-hdf5, but this is NOT recommended: all current output
formats are built on HDF5, and the legacy ASCII output is unmaintained.])])

  AC_MSG_CHECKING([for HDF5 Fortran compile and link flags])
  hdf5_show=`"$H5FC" -show 2>/dev/null`
  AS_IF([test -z "$hdf5_show"],
    [AC_MSG_FAILURE([$H5FC -show returned no compiler or linker flags])])

  hdf5_wrapper_flags=${hdf5_show#* }
  for hdf5_flag in $hdf5_wrapper_flags
  do
    AS_CASE([$hdf5_flag],
      [-show], [],
      [-I*|-D*|-U*], [HDF5_FCFLAGS="$HDF5_FCFLAGS $hdf5_flag"],
      [HDF5_LIBS="$HDF5_LIBS $hdf5_flag"])
  done
  AS_IF([test -z "$HDF5_FCFLAGS" || test -z "$HDF5_LIBS"],
    [AC_MSG_FAILURE([could not parse compile and link flags from: $hdf5_show])])
  AC_MSG_RESULT([yes])

  hdf5_save_FCFLAGS=$FCFLAGS
  hdf5_save_LIBS=$LIBS
  FCFLAGS="$FCFLAGS $HDF5_FCFLAGS"
  LIBS="$LIBS $HDF5_LIBS"
  AC_MSG_CHECKING([whether $FC can compile and link the HDF5 Fortran module])
  AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
      program conftest
      use hdf5
      implicit none
      integer :: hdferr
      call h5open_f(hdferr)
      call h5close_f(hdferr)
      end program conftest
    ]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_FAILURE(
       [the HDF5 Fortran flags reported by $H5FC do not work with $FC.
Ensure h5fc and the simulation Fortran compiler use compatible ABIs.
h5fc -show: $hdf5_show])])
  FCFLAGS=$hdf5_save_FCFLAGS
  LIBS=$hdf5_save_LIBS

  hdf5_showconfig=`LC_ALL=C "$H5FC" -showconfig 2>/dev/null`
  hdf5_parallel_value=`printf '%s\n' "$hdf5_showconfig" |
    sed -n 's/^[[[:space:]]]*Parallel HDF5:[[[:space:]]]*//p' |
    sed -n '1{s/[[[:space:]]]*$//;p;}' |
    LC_ALL=C tr '[[:upper:]]' '[[:lower:]]'`
  AS_CASE([$hdf5_parallel_value],
    [yes|on], [HDF5_PARALLEL=yes],
    [no|off], [HDF5_PARALLEL=no],
    [HDF5_PARALLEL=unknown])

  AS_IF([test "x$HDF5_PARALLEL" = xyes],
    [AC_MSG_WARN(
       [parallel HDF5 detected through $H5FC. NBODY6++GPU HDF5 output supports serial HDF5 only.
Load a serial HDF5 module or put its h5fc first in PATH.])],
    [test "x$HDF5_PARALLEL" = xunknown],
    [AC_MSG_WARN(
       [could not confirm that the HDF5 library discovered through $H5FC is serial.
NBODY6++GPU HDF5 output supports serial HDF5 only.])])
])


# FIXME: Replace `main' with a function in `-lmpi':
# AC_CHECK_LIB([mpi], [main])


#C++ part
AC_LANG([C++])
local_cxxflags=$CXXFLAGS
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_LIB([stdc++],[main],[LIBS=$LIBS' -lstdc++'])

#AS_IF([test "x$enable_debug" != xno],
#      [CXXFLAGS=$local_cxxflags" -O2"],
#      [CXXFLAGS=$local_cxxflags" -O3"])
CXXFLAGS=$local_cxxflags" -O3"

AC_OPENMP

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

AS_IF([test "x$enable_mcmodel" != xno],
   [AX_CHECK_COMPILE_FLAG([-mcmodel=$enable_mcmodel],
      [CXXFLAGS=$CXXFLAGS" -mcmodel="$enable_mcmodel],
      [AC_MSG_FAILURE(
          [$CXX compiler option -mcmodel=$enable_mcmodel is not avaiable])])])

# SIMD part
# Handle SIMDe for ARM/cross-platform SIMD support
AS_IF([test "x$with_simde" != xno],
   [AC_MSG_CHECKING([for SIMDe support])
    # Check if simde path was provided
    AS_IF([test "x$with_simde" != xyes && test "x$with_simde" != x],
       [AS_IF([test ! -d "$with_simde"],
          [AC_MSG_FAILURE([SIMDe include location does not exist: $with_simde
Install SIMDe v0.8.2 with:
  git clone --branch v0.8.2 --depth 1 https://github.com/simd-everywhere/simde.git "$simde_default_path"
ARM64 example: $arm_configure_example])])
        for simde_header in sse.h sse2.h sse3.h sse4.1.h; do
          AS_IF([test ! -f "$with_simde/simde/x86/$simde_header"],
            [AC_MSG_FAILURE([SIMDe header is missing: $with_simde/simde/x86/$simde_header
Expected include location: $with_simde
Install SIMDe v0.8.2 with:
  git clone --branch v0.8.2 --depth 1 https://github.com/simd-everywhere/simde.git "$simde_default_path"
ARM64 example: $arm_configure_example])])
        done
        SIMDE_CXXFLAGS="-I$with_simde"
        AC_MSG_RESULT([using path: $with_simde])],
       [SIMDE_CXXFLAGS=""
        AS_IF([test "x$with_simde" = x], [with_simde=yes])
        AC_MSG_RESULT([using system include paths])])

    # Add SIMDe flags
    CXXFLAGS=$CXXFLAGS" -DNBODY_USE_SIMDE $SIMDE_CXXFLAGS"

    AC_MSG_CHECKING([whether the required SIMDe headers compile])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <simde/x86/sse.h>
#include <simde/x86/sse2.h>
#include <simde/x86/sse3.h>
#include <simde/x86/sse4.1.h>
]], [[simde__m128 value = simde_mm_setzero_ps(); (void) value;]])],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_FAILURE([the required SIMDe headers are not compilable.
Checked include location: ${with_simde:-system compiler include paths}
Required headers: simde/x86/sse.h, sse2.h, sse3.h, sse4.1.h
Install SIMDe v0.8.2 with:
  git clone --branch v0.8.2 --depth 1 https://github.com/simd-everywhere/simde.git "$simde_default_path"
ARM64 example: $arm_configure_example])])

    # Allow SSE on ARM with SIMDe
    AS_CASE([$target_architecture],
       [ARM64], [
          SIMD_CHECK=1])
   ])

AS_IF([test "x$OPENMP_CXXFLAGS" != x],
   [AS_IF([test "x$enable_simd" = xtest],
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
               [# Skip -msse/-mavx flag check for SIMDe on ARM
                AS_IF([test "x$with_simde" != xno],
                   [AS_CASE([$target_architecture],
                      [ARM64], [
                         # SIMDe on ARM - don't need -msse flag, RESULT set here
                         AC_MSG_NOTICE([Using SIMDe for SIMD support on ARM])
                         RESULT=$RESULT".sse"
                         SIMD_CHECK=1],
                      [# SIMDe on x86 - still check for SSE support
                       AC_CHECK_FILE(["/proc/cpuinfo"],
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
                             [$CXX compiler option -m$enable_simd is not available])])])],
                   [# No SIMDe - regular x86 SIMD
                    AC_CHECK_FILE(["/proc/cpuinfo"],
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
                          [$CXX compiler option -m$enable_simd is not available])])])])])],
     [AS_IF([test "x$with_simde" = xno],
        [enable_simd=no
         SIMD_CHECK=-2])])

AS_IF([test "x$enable_simd" = xno],
   [AC_MSG_WARN([SIMD is disabled (--enable-simd=no); OpenMP is being disabled
automatically for this build too. The plain Fortran nbint.F force routine
used without SIMD is not thread-safe under OpenMP (see README Known
Problems #3), so this build will be single-threaded. Use only for
debugging.])
    OPENMP_FCFLAGS=
    OPENMP_CXXFLAGS=])

# FIXME: Replace `main' with a function in `-llvisit':
# AC_CHECK_LIB([lvisit], [main])
# FIXME: Replace `main' with a function in `-llvisit_nbody':
# AC_CHECK_LIB([lvisit_nbody], [main])
# FIXME: Replace `main' with a function in `-lvisit':
# AC_CHECK_LIB([visit], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday memset sqrt strchr])

# CUDA/HIP GPU backends
AC_ARG_VAR([NVCC], [CUDA compiler command])
AC_ARG_VAR([HIPCC], [HIP compiler command])
AC_ARG_WITH([cuda],
     [AS_HELP_STRING([--with-cuda=PREFIX],
         [Prefix of your CUDA installation])],
     [cuda_prefix=$withval],
     [cuda_prefix=""])

AC_ARG_WITH([hip],
     [AS_HELP_STRING([--with-hip=PREFIX],
         [Prefix of a standard AMD ROCm installation])],
     [hip_prefix=$withval],
     [hip_prefix=""])

AC_ARG_WITH([cuda_sdk],
     [AS_HELP_STRING([--with-cuda-sdk=PREFIX],
         [Prefix of your CUDA SDKinstallation])],
     [cuda_sdk_prefix=$withval],
     [cuda_sdk_prefix="/usr/local/cuda/samples"])

# Normalize relative CUDA paths to absolute paths (build runs in ./build)
AS_CASE([$cuda_prefix],
   [/*], [],
   [""], [],
   [*], [AS_IF([test -d "$cuda_prefix"],
          [cuda_prefix=`cd "$cuda_prefix" && pwd`])])
AS_CASE([$hip_prefix],
   [/*], [],
   [""], [],
   [*], [AS_IF([test -d "$hip_prefix"],
          [hip_prefix=`cd "$hip_prefix" && pwd`])])
AS_CASE([$cuda_sdk_prefix],
   [/*], [],
   [*], [AS_IF([test -d "$cuda_sdk_prefix"],
          [cuda_sdk_prefix=`cd "$cuda_sdk_prefix" && pwd`])])

GPU_BACKEND=none
GPU_COMPILER=
GPU_PREFIX=
GPUFLAGS=""
AS_IF([test "x$enable_gpu" != xno], [
  AS_IF([test -z "$NVCC"], [
    AS_IF([test -n "$cuda_prefix" && test -x "$cuda_prefix/bin/nvcc"], [NVCC="$cuda_prefix/bin/nvcc"],
      [AC_PATH_PROG([NVCC], [nvcc], [])])])
  AS_IF([test -z "$HIPCC"], [
    AS_IF([test -n "$hip_prefix" && test -x "$hip_prefix/bin/hipcc"], [HIPCC="$hip_prefix/bin/hipcc"], [
      AC_PATH_PROG([HIPCC], [hipcc], [], [/opt/rocm/bin:/usr/local/rocm/bin:$PATH])])])
  AS_CASE([$gpu_backend],
    [cuda], [test -n "$NVCC" || AC_MSG_FAILURE([CUDA selected but nvcc was not found; set NVCC or use --with-cuda=PREFIX]); GPU_BACKEND=cuda],
    [hip],  [test -n "$HIPCC" || AC_MSG_FAILURE([HIP selected but hipcc was not found; set HIPCC or use --with-hip=PREFIX]); GPU_BACKEND=hip],
    [auto], [AS_IF([test -n "$NVCC"], [GPU_BACKEND=cuda], [test -n "$HIPCC"], [GPU_BACKEND=hip],
                    [AC_MSG_FAILURE([No CUDA or HIP compiler found; install nvcc/hipcc or configure with --disable-gpu])])])

  AS_IF([test "x$GPU_BACKEND" = xcuda], [
    GPU_COMPILER=$NVCC
    GPUFLAGS="-O3 -Xcompiler=-fPIC"
    AS_IF([test "x$OPENMP_CXXFLAGS" != x], [GPUFLAGS="$GPUFLAGS -Xcompiler=$OPENMP_CXXFLAGS"])
    AS_IF([test -z "$cuda_prefix"], [cuda_bindir=`dirname "$NVCC"`; cuda_prefix=`cd "$cuda_bindir/.." 2>/dev/null && pwd`])
    GPU_PREFIX=$cuda_prefix
    GPU_LINK_FLAGS="-L$GPU_PREFIX/lib64 -lcudart"
    AC_MSG_CHECKING([whether CUDA compiler and runtime link])
    cat > conftest.cu <<'_ACEOF'
#include <cuda_runtime.h>
int main(void) { return cudaFree(0) == cudaSuccess ? 0 : 1; }
_ACEOF
    AS_IF(["$NVCC" $GPUFLAGS conftest.cu $LDFLAGS $GPU_LINK_FLAGS -o conftest$ac_exeext >/dev/null 2>&1],
      [AC_MSG_RESULT([yes])], [AC_MSG_RESULT([no]); rm -f conftest*; AC_MSG_FAILURE([CUDA compile/link test failed for $NVCC and $GPU_PREFIX])])
    rm -f conftest*
    CUFLAGS="$GPUFLAGS -D CUDA_5 -I ../extra_inc/cuda -I ../include"
    AS_IF([test "x$OPENMP_CXXFLAGS" != x],
      [EXTRAOBJ=$EXTRAOBJ' $(CUDA_OBJECTS)'],
      [EXTRAOBJ=$EXTRAOBJ' $(CUDA_OBJECTS_SINGLE)'])
    RESULT=$RESULT.gpu
  ], [
    GPU_COMPILER=$HIPCC
    GPUFLAGS="-O3 -fPIC"
dnl OPENMP_CXXFLAGS was probed with the host C++ compiler (g++ -> libgomp), but
dnl hipcc is clang-based and resolves a bare -fopenmp to libomp.  The final link
dnl is done by the Fortran compiler, which only pulls in libgomp, so a bare
dnl -fopenmp here leaves __kmpc_fork_call and friends undefined.  Pin clang to
dnl the GNU OpenMP ABI so both halves of the program share one runtime.
    AS_IF([test "x$OPENMP_CXXFLAGS" = "x-fopenmp"],
          [hip_openmp_flags="-fopenmp=libgomp"],
          [hip_openmp_flags=$OPENMP_CXXFLAGS])
    AS_IF([test "x$hip_openmp_flags" != x], [GPUFLAGS="$GPUFLAGS $hip_openmp_flags"])
    AS_IF([test -z "$hip_prefix"], [hip_bindir=`dirname "$HIPCC"`; hip_prefix=`cd "$hip_bindir/.." 2>/dev/null && pwd`])
    GPU_PREFIX=$hip_prefix
    GPU_LINK_FLAGS="-L$GPU_PREFIX/lib -L$GPU_PREFIX/lib64 -lamdhip64"
    AC_MSG_CHECKING([whether HIP compiler and runtime link])
    cat > conftest.hip <<'_ACEOF'
#include <hip/hip_runtime.h>
int main(void) { return hipFree(0) == hipSuccess ? 0 : 1; }
_ACEOF
    AS_IF(["$HIPCC" $GPUFLAGS conftest.hip $LDFLAGS $GPU_LINK_FLAGS -o conftest$ac_exeext >/dev/null 2>&1],
      [AC_MSG_RESULT([yes])], [AC_MSG_RESULT([no]); rm -f conftest*; AC_MSG_FAILURE([HIP compile/link test failed for $HIPCC and $GPU_PREFIX])])
    rm -f conftest*
    HIPFLAGS="$GPUFLAGS -I../src/Main"
    AS_IF([test "x$OPENMP_CXXFLAGS" != x],
      [EXTRAOBJ=$EXTRAOBJ' $(HIP_OBJECTS)'],
      [EXTRAOBJ=$EXTRAOBJ' $(HIP_OBJECTS_SINGLE)'])
    RESULT=$RESULT.hip
  ])
  LIBS="$GPU_LINK_FLAGS $LIBS"
  FCFLAGS=$FCFLAGS' $(GPU_FLAGS)'
  EXTRASRC=$EXTRASRC' $(GPU_FSOURCES)'
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
      FCFLAGS=$FCFLAGS' '$OPENMP_FCFLAGS' ${OMP_FLAGS}',
      FCFLAGS=$FCFLAGS' -I../extra_inc/noomp')

AS_IF([test "x$OPENMP_CXXFLAGS" != x ],
      CXXFLAGS=$CXXFLAGS' '$OPENMP_CXXFLAGS' ${OMP_FLAGS}',
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
         [Size of basic parameters (NMAX, KMAX, LMAX, MMAX) in params.h. Available values: 128m, 64m, 32m, 16m, 1m, 512k, 256k, 128k, 64k, 32k, 16k, 8k, 4k, 2k, 1k, custom. The value of this size is NMAX, with prefix 'b' the K.S. binary size KMAX is half of NMAX and NMAX is set to 1.5*NMAX. LMAX is 600 above 64k, 500 above 16k. MMAX is 1024 above 64k, 512 above 16k. MMAX is doubled with prefix 'b'. Defaultsize is b16k])],
     [par_size=$withval],
     [par_size=b16k])

AS_IF([test "x$par_size" == x128m],
      [NMAX=134217728
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb128m],
      [NMAX=201326592
       KMAX=67108864
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x64m],
      [NMAX=67108864
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb64m],
      [NMAX=100663296
       KMAX=33554432
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x32m],
      [NMAX=33554432
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb32m],
      [NMAX=50331648
       KMAX=16777216
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x16m],
      [NMAX=16777216
       KMAX=2097152
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
      [NMAX=15107200
       KMAX=5120000
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x8m],
      [NMAX=8388608
       KMAX=131072
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb8m],
      [NMAX=6291456
       KMAX=2097152
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x4m],
      [NMAX=4194304
       KMAX=524288
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
      [test "x$par_size" == x4k],
      [NMAX=4096
       KMAX=1024
       LMAX=400
       MMAX=128],
      [test "x$par_size" == xb4k],
      [NMAX=6144
       KMAX=2048
       LMAX=400
       MMAX=256],
      [test "x$par_size" == x2k],
      [NMAX=2048
       KMAX=512
       LMAX=400
       MMAX=128],
      [test "x$par_size" == xb2k],
      [NMAX=3072
       KMAX=1024
       LMAX=400
       MMAX=256],
      [test "x$par_size" == x1k],
      [NMAX=1024
       KMAX=256
       LMAX=400
       MMAX=128],
      [test "x$par_size" == xb1k],
      [NMAX=1536
       KMAX=512
       LMAX=400
       MMAX=256],
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


# Extral tools
AS_IF([test "x$enable_tools" != xno],
   [EXTRARESULT=' $(EXTRATOOLS)'],
   [EXTRARESULT=''])

AC_SUBST([CUFLAGS])
AC_SUBST([HIPFLAGS])
AC_SUBST([NVCC])
AC_SUBST([HIPCC])
AC_SUBST([FCFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([RESULT])
AC_SUBST([EXTRASRC])
AC_SUBST([EXTRAOBJ])
AC_SUBST([EXTRARESULT])
AC_SUBST([MPI_FLAGS])
AC_SUBST([HDF5_FCFLAGS])
AC_SUBST([HDF5_LIBS])
AC_SUBST([HDF5_PARALLEL])

AC_SUBST([NMAX])
AC_SUBST([KMAX])
AC_SUBST([LMAX])
AC_SUBST([MMAX])

AC_MSG_NOTICE([----------------------Configure Summary--------------------])
AC_MSG_NOTICE([     Target architecture: $target_architecture])
AC_MSG_NOTICE([     Canonical host:      $host])
AS_IF([test "x$build" != "x$host"], [
  AC_MSG_NOTICE([     Canonical build:     $build])
])
AC_MSG_NOTICE([--Features:])
AC_MSG_NOTICE([     Using MPI:      $enable_mpi])
AC_MSG_NOTICE([     Using GPU:      $enable_gpu])
AC_MSG_NOTICE([     GPU backend:    $GPU_BACKEND])
AC_MSG_NOTICE([     GPU prefix:     $GPU_PREFIX])
AC_MSG_NOTICE([     GPU flags:      $GPUFLAGS])
AC_MSG_NOTICE([     Using SIMD:     $enable_simd])
AS_IF([test "x$with_simde" != xno],
   [AS_IF([test "x$with_simde" = xyes],
      [AC_MSG_NOTICE([     Using SIMDe:    yes (system path)])],
      [AC_MSG_NOTICE([     Using SIMDe:    yes (path: $with_simde)])])],
   [AC_MSG_NOTICE([     Using SIMDe:    no])])
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
AS_IF([test "x$enable_hdf5" != xno], [
  AC_MSG_NOTICE([     HDF5 wrapper:   $H5FC])
  AC_MSG_NOTICE([     Parallel HDF5:  $HDF5_PARALLEL])
])
AC_MSG_NOTICE([     TT package:     $enable_tt])
AC_MSG_NOTICE([--Compilers:])
AC_MSG_NOTICE([     C++ compiler:      $CXX])
AC_MSG_NOTICE([     Fortran compiler:  $FC])
AC_MSG_NOTICE([     GPU compiler:      $GPU_COMPILER])
AC_MSG_NOTICE([--Targets:])
AC_MSG_NOTICE([     Program name:   nbody6++$RESULT])
AC_MSG_NOTICE([     Extra tools:    $enable_tools])
AC_MSG_NOTICE([     Install path:   $prefix])
AC_MSG_NOTICE([--Parameters:])
AC_MSG_NOTICE([     NMAX:           $NMAX])
AC_MSG_NOTICE([     KMAX:           $KMAX])
AC_MSG_NOTICE([     LMAX:           $LMAX])
AC_MSG_NOTICE([     MMAX:           $MMAX])
AC_MSG_NOTICE([-----------------------------------------------------------])
AS_IF([test "x$enable_tools" != xno],
   [AC_MSG_NOTICE([  Please add $prefix/lib into your runtime load library path:])
    AC_MSG_NOTICE([     LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$prefix/lib])
    AC_MSG_NOTICE([  in your pre-configure script of shell (such as ~/.bashrc if you use bash)])
    AC_MSG_NOTICE([-----------------------------------------------------------])])
AC_MSG_NOTICE([Reference paper: Wang, L., et al., 2015, MNRAS, 450, 4070 (http://arxiv.org/abs/1504.03687)])
AC_MSG_NOTICE([Email: Long Wang<longwang.astro@gmail.com>])
AC_MSG_NOTICE([-----------------------------------------------------------])
AC_CONFIG_FILES([Makefile build/Makefile include/params.h])
AC_OUTPUT
