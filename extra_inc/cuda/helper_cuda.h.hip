/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

////////////////////////////////////////////////////////////////////////////////
// These are CUDA Helper functions for initialization and error checking

#ifndef HELPER_CUDA_H
#define HELPER_CUDA_H

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <helper_string.h.hip>

//#include <string>
//#include <iostream>
//#include <sstream>

// Note, it is required that your SDK sample to include the proper header files, please
// refer the CUDA examples for examples of the needed CUDA headers, which may change depending
// on which CUDA functions are used.

// CUDA Runtime error messages
#ifdef __DRIVER_TYPES_H__
static const char *_cudaGetErrorEnum(hipError_t error)
{
    switch (error)
    {
        case hipSuccess:
            return "hipSuccess";

        case hipErrorMissingConfiguration:
            return "hipErrorMissingConfiguration";

        case hipErrorOutOfMemory:
            return "hipErrorOutOfMemory";

        case hipErrorNotInitialized:
            return "hipErrorNotInitialized";

        case hipErrorLaunchFailure:
            return "hipErrorLaunchFailure";

        case hipErrorPriorLaunchFailure:
            return "hipErrorPriorLaunchFailure";

        case hipErrorLaunchTimeOut:
            return "hipErrorLaunchTimeOut";

        case hipErrorLaunchOutOfResources:
            return "hipErrorLaunchOutOfResources";

        case hipErrorInvalidDeviceFunction:
            return "hipErrorInvalidDeviceFunction";

        case hipErrorInvalidConfiguration:
            return "hipErrorInvalidConfiguration";

        case hipErrorInvalidDevice:
            return "hipErrorInvalidDevice";

        case hipErrorInvalidValue:
            return "hipErrorInvalidValue";

        case hipErrorInvalidPitchValue:
            return "hipErrorInvalidPitchValue";

        case hipErrorInvalidSymbol:
            return "hipErrorInvalidSymbol";

        case hipErrorMapFailed:
            return "hipErrorMapFailed";

        case hipErrorUnmapFailed:
            return "hipErrorUnmapFailed";

        case cudaErrorInvalidHostPointer:
            return "cudaErrorInvalidHostPointer";

        case hipErrorInvalidDevicePointer:
            return "hipErrorInvalidDevicePointer";

        case hipErrorInvalidTexture:
            return "hipErrorInvalidTexture";

        case cudaErrorInvalidTextureBinding:
            return "cudaErrorInvalidTextureBinding";

        case hipErrorInvalidChannelDescriptor:
            return "hipErrorInvalidChannelDescriptor";

        case hipErrorInvalidMemcpyDirection:
            return "hipErrorInvalidMemcpyDirection";

        case cudaErrorAddressOfConstant:
            return "cudaErrorAddressOfConstant";

        case cudaErrorTextureFetchFailed:
            return "cudaErrorTextureFetchFailed";

        case cudaErrorTextureNotBound:
            return "cudaErrorTextureNotBound";

        case cudaErrorSynchronizationError:
            return "cudaErrorSynchronizationError";

        case cudaErrorInvalidFilterSetting:
            return "cudaErrorInvalidFilterSetting";

        case cudaErrorInvalidNormSetting:
            return "cudaErrorInvalidNormSetting";

        case cudaErrorMixedDeviceExecution:
            return "cudaErrorMixedDeviceExecution";

        case hipErrorDeinitialized:
            return "hipErrorDeinitialized";

        case hipErrorUnknown:
            return "hipErrorUnknown";

        case cudaErrorNotYetImplemented:
            return "cudaErrorNotYetImplemented";

        case cudaErrorMemoryValueTooLarge:
            return "cudaErrorMemoryValueTooLarge";

        case hipErrorInvalidHandle:
            return "hipErrorInvalidHandle";

        case hipErrorNotReady:
            return "hipErrorNotReady";

        case hipErrorInsufficientDriver:
            return "hipErrorInsufficientDriver";

        case hipErrorSetOnActiveProcess:
            return "hipErrorSetOnActiveProcess";

        case cudaErrorInvalidSurface:
            return "cudaErrorInvalidSurface";

        case hipErrorNoDevice:
            return "hipErrorNoDevice";

        case hipErrorECCNotCorrectable:
            return "hipErrorECCNotCorrectable";

        case hipErrorSharedObjectSymbolNotFound:
            return "hipErrorSharedObjectSymbolNotFound";

        case hipErrorSharedObjectInitFailed:
            return "hipErrorSharedObjectInitFailed";

        case hipErrorUnsupportedLimit:
            return "hipErrorUnsupportedLimit";

        case cudaErrorDuplicateVariableName:
            return "cudaErrorDuplicateVariableName";

        case cudaErrorDuplicateTextureName:
            return "cudaErrorDuplicateTextureName";

        case cudaErrorDuplicateSurfaceName:
            return "cudaErrorDuplicateSurfaceName";

        case cudaErrorDevicesUnavailable:
            return "cudaErrorDevicesUnavailable";

        case hipErrorInvalidImage:
            return "hipErrorInvalidImage";

        case hipErrorNoBinaryForGpu:
            return "hipErrorNoBinaryForGpu";

        case cudaErrorIncompatibleDriverContext:
            return "cudaErrorIncompatibleDriverContext";

        case hipErrorPeerAccessAlreadyEnabled:
            return "hipErrorPeerAccessAlreadyEnabled";

        case hipErrorPeerAccessNotEnabled:
            return "hipErrorPeerAccessNotEnabled";

        case hipErrorContextAlreadyInUse:
            return "hipErrorContextAlreadyInUse";

        case hipErrorProfilerDisabled:
            return "hipErrorProfilerDisabled";

        case hipErrorProfilerNotInitialized:
            return "hipErrorProfilerNotInitialized";

        case hipErrorProfilerAlreadyStarted:
            return "hipErrorProfilerAlreadyStarted";

        case hipErrorProfilerAlreadyStopped:
            return "hipErrorProfilerAlreadyStopped";

#if __CUDA_API_VERSION >= 0x4000

        case hipErrorAssert:
            return "hipErrorAssert";

        case cudaErrorTooManyPeers:
            return "cudaErrorTooManyPeers";

        case hipErrorHostMemoryAlreadyRegistered:
            return "hipErrorHostMemoryAlreadyRegistered";

        case hipErrorHostMemoryNotRegistered:
            return "hipErrorHostMemoryNotRegistered";
#endif

        case cudaErrorStartupFailure:
            return "cudaErrorStartupFailure";

        case cudaErrorApiFailureBase:
            return "cudaErrorApiFailureBase";
    }

    return "<unknown>";
}
#endif

#ifdef __cuda_cuda_h__
// CUDA Driver API errors
static const char *_cudaGetErrorEnum(hipError_t error)
{
    switch (error)
    {
        case hipSuccess:
            return "hipSuccess";

        case hipErrorInvalidValue:
            return "hipErrorInvalidValue";

        case hipErrorOutOfMemory:
            return "hipErrorOutOfMemory";

        case hipErrorNotInitialized:
            return "hipErrorNotInitialized";

        case hipErrorDeinitialized:
            return "hipErrorDeinitialized";

        case hipErrorProfilerDisabled:
            return "hipErrorProfilerDisabled";

        case hipErrorProfilerNotInitialized:
            return "hipErrorProfilerNotInitialized";

        case hipErrorProfilerAlreadyStarted:
            return "hipErrorProfilerAlreadyStarted";

        case hipErrorProfilerAlreadyStopped:
            return "hipErrorProfilerAlreadyStopped";

        case hipErrorNoDevice:
            return "hipErrorNoDevice";

        case hipErrorInvalidDevice:
            return "hipErrorInvalidDevice";

        case hipErrorInvalidImage:
            return "hipErrorInvalidImage";

        case hipErrorInvalidContext:
            return "hipErrorInvalidContext";

        case hipErrorContextAlreadyCurrent:
            return "hipErrorContextAlreadyCurrent";

        case hipErrorMapFailed:
            return "hipErrorMapFailed";

        case hipErrorUnmapFailed:
            return "hipErrorUnmapFailed";

        case hipErrorArrayIsMapped:
            return "hipErrorArrayIsMapped";

        case hipErrorAlreadyMapped:
            return "hipErrorAlreadyMapped";

        case hipErrorNoBinaryForGpu:
            return "hipErrorNoBinaryForGpu";

        case hipErrorAlreadyAcquired:
            return "hipErrorAlreadyAcquired";

        case hipErrorNotMapped:
            return "hipErrorNotMapped";

        case hipErrorNotMappedAsArray:
            return "hipErrorNotMappedAsArray";

        case hipErrorNotMappedAsPointer:
            return "hipErrorNotMappedAsPointer";

        case hipErrorECCNotCorrectable:
            return "hipErrorECCNotCorrectable";

        case hipErrorUnsupportedLimit:
            return "hipErrorUnsupportedLimit";

        case hipErrorContextAlreadyInUse:
            return "hipErrorContextAlreadyInUse";

        case hipErrorInvalidSource:
            return "hipErrorInvalidSource";

        case hipErrorFileNotFound:
            return "hipErrorFileNotFound";

        case hipErrorSharedObjectSymbolNotFound:
            return "hipErrorSharedObjectSymbolNotFound";

        case hipErrorSharedObjectInitFailed:
            return "hipErrorSharedObjectInitFailed";

        case hipErrorOperatingSystem:
            return "hipErrorOperatingSystem";

        case hipErrorInvalidHandle:
            return "hipErrorInvalidHandle";

        case hipErrorNotFound:
            return "hipErrorNotFound";

        case hipErrorNotReady:
            return "hipErrorNotReady";

        case hipErrorLaunchFailure:
            return "hipErrorLaunchFailure";

        case hipErrorLaunchOutOfResources:
            return "hipErrorLaunchOutOfResources";

        case hipErrorLaunchTimeOut:
            return "hipErrorLaunchTimeOut";

        case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING:
            return "CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING";

        case hipErrorPeerAccessAlreadyEnabled:
            return "hipErrorPeerAccessAlreadyEnabled";

        case hipErrorPeerAccessNotEnabled:
            return "hipErrorPeerAccessNotEnabled";

        case hipErrorSetOnActiveProcess:
            return "hipErrorSetOnActiveProcess";

        case hipErrorContextIsDestroyed:
            return "hipErrorContextIsDestroyed";

        case hipErrorAssert:
            return "hipErrorAssert";

        case CUDA_ERROR_TOO_MANY_PEERS:
            return "CUDA_ERROR_TOO_MANY_PEERS";

        case hipErrorHostMemoryAlreadyRegistered:
            return "hipErrorHostMemoryAlreadyRegistered";

        case hipErrorHostMemoryNotRegistered:
            return "hipErrorHostMemoryNotRegistered";

        case hipErrorUnknown:
            return "hipErrorUnknown";
    }

    return "<unknown>";
}
#endif

#ifdef CUBLAS_API_H_
// cuBLAS API errors
static const char *_cudaGetErrorEnum(hipblasStatus_t error)
{
    switch (error)
    {
        case HIPBLAS_STATUS_SUCCESS:
            return "HIPBLAS_STATUS_SUCCESS";

        case HIPBLAS_STATUS_NOT_INITIALIZED:
            return "HIPBLAS_STATUS_NOT_INITIALIZED";

        case HIPBLAS_STATUS_ALLOC_FAILED:
            return "HIPBLAS_STATUS_ALLOC_FAILED";

        case HIPBLAS_STATUS_INVALID_VALUE:
            return "HIPBLAS_STATUS_INVALID_VALUE";

        case HIPBLAS_STATUS_ARCH_MISMATCH:
            return "HIPBLAS_STATUS_ARCH_MISMATCH";

        case HIPBLAS_STATUS_MAPPING_ERROR:
            return "HIPBLAS_STATUS_MAPPING_ERROR";

        case HIPBLAS_STATUS_EXECUTION_FAILED:
            return "HIPBLAS_STATUS_EXECUTION_FAILED";

        case HIPBLAS_STATUS_INTERNAL_ERROR:
            return "HIPBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}
#endif

#ifdef _CUFFT_H_
// cuFFT API errors
static const char *_cudaGetErrorEnum(hipfftResult error)
{
    switch (error)
    {
        case HIPFFT_SUCCESS:
            return "HIPFFT_SUCCESS";

        case HIPFFT_INVALID_PLAN:
            return "HIPFFT_INVALID_PLAN";

        case HIPFFT_ALLOC_FAILED:
            return "HIPFFT_ALLOC_FAILED";

        case HIPFFT_INVALID_TYPE:
            return "HIPFFT_INVALID_TYPE";

        case HIPFFT_INVALID_VALUE:
            return "HIPFFT_INVALID_VALUE";

        case HIPFFT_INTERNAL_ERROR:
            return "HIPFFT_INTERNAL_ERROR";

        case HIPFFT_EXEC_FAILED:
            return "HIPFFT_EXEC_FAILED";

        case HIPFFT_SETUP_FAILED:
            return "HIPFFT_SETUP_FAILED";

        case HIPFFT_INVALID_SIZE:
            return "HIPFFT_INVALID_SIZE";

        case HIPFFT_UNALIGNED_DATA:
            return "HIPFFT_UNALIGNED_DATA";
    }

    return "<unknown>";
}
#endif


#ifdef CUSPARSEAPI
// cuSPARSE API errors
static const char *_cudaGetErrorEnum(hipsparseStatus_t error)
{
    switch (error)
    {
        case HIPSPARSE_STATUS_SUCCESS:
            return "HIPSPARSE_STATUS_SUCCESS";

        case HIPSPARSE_STATUS_NOT_INITIALIZED:
            return "HIPSPARSE_STATUS_NOT_INITIALIZED";

        case HIPSPARSE_STATUS_ALLOC_FAILED:
            return "HIPSPARSE_STATUS_ALLOC_FAILED";

        case HIPSPARSE_STATUS_INVALID_VALUE:
            return "HIPSPARSE_STATUS_INVALID_VALUE";

        case HIPSPARSE_STATUS_ARCH_MISMATCH:
            return "HIPSPARSE_STATUS_ARCH_MISMATCH";

        case HIPSPARSE_STATUS_MAPPING_ERROR:
            return "HIPSPARSE_STATUS_MAPPING_ERROR";

        case HIPSPARSE_STATUS_EXECUTION_FAILED:
            return "HIPSPARSE_STATUS_EXECUTION_FAILED";

        case HIPSPARSE_STATUS_INTERNAL_ERROR:
            return "HIPSPARSE_STATUS_INTERNAL_ERROR";

        case HIPSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            return "HIPSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    }

    return "<unknown>";
}
#endif

#ifdef CURAND_H_
// cuRAND API errors
static const char *_cudaGetErrorEnum(hiprandStatus_t error)
{
    switch (error)
    {
        case HIPRAND_STATUS_SUCCESS:
            return "HIPRAND_STATUS_SUCCESS";

        case HIPRAND_STATUS_VERSION_MISMATCH:
            return "HIPRAND_STATUS_VERSION_MISMATCH";

        case HIPRAND_STATUS_NOT_INITIALIZED:
            return "HIPRAND_STATUS_NOT_INITIALIZED";

        case HIPRAND_STATUS_ALLOCATION_FAILED:
            return "HIPRAND_STATUS_ALLOCATION_FAILED";

        case HIPRAND_STATUS_TYPE_ERROR:
            return "HIPRAND_STATUS_TYPE_ERROR";

        case HIPRAND_STATUS_OUT_OF_RANGE:
            return "HIPRAND_STATUS_OUT_OF_RANGE";

        case HIPRAND_STATUS_LENGTH_NOT_MULTIPLE:
            return "HIPRAND_STATUS_LENGTH_NOT_MULTIPLE";

        case HIPRAND_STATUS_DOUBLE_PRECISION_REQUIRED:
            return "HIPRAND_STATUS_DOUBLE_PRECISION_REQUIRED";

        case HIPRAND_STATUS_LAUNCH_FAILURE:
            return "HIPRAND_STATUS_LAUNCH_FAILURE";

        case HIPRAND_STATUS_PREEXISTING_FAILURE:
            return "HIPRAND_STATUS_PREEXISTING_FAILURE";

        case HIPRAND_STATUS_INITIALIZATION_FAILED:
            return "HIPRAND_STATUS_INITIALIZATION_FAILED";

        case HIPRAND_STATUS_ARCH_MISMATCH:
            return "HIPRAND_STATUS_ARCH_MISMATCH";

        case HIPRAND_STATUS_INTERNAL_ERROR:
            return "HIPRAND_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}
#endif

#ifdef NV_NPPIDEFS_H
// NPP API errors
static const char *_cudaGetErrorEnum(NppStatus error)
{
    switch (error)
    {
        case NPP_NOT_SUPPORTED_MODE_ERROR:
            return "NPP_NOT_SUPPORTED_MODE_ERROR";

        case NPP_ROUND_MODE_NOT_SUPPORTED_ERROR:
            return "NPP_ROUND_MODE_NOT_SUPPORTED_ERROR";

        case NPP_RESIZE_NO_OPERATION_ERROR:
            return "NPP_RESIZE_NO_OPERATION_ERROR";

        case NPP_NOT_SUFFICIENT_COMPUTE_CAPABILITY:
            return "NPP_NOT_SUFFICIENT_COMPUTE_CAPABILITY";

        case NPP_BAD_ARG_ERROR:
            return "NPP_BAD_ARG_ERROR";

        case NPP_LUT_NUMBER_OF_LEVELS_ERROR:
            return "NPP_LUT_NUMBER_OF_LEVELS_ERROR";

        case NPP_TEXTURE_BIND_ERROR:
            return "NPP_TEXTURE_BIND_ERROR";

        case NPP_COEFF_ERROR:
            return "NPP_COEFF_ERROR";

        case NPP_RECT_ERROR:
            return "NPP_RECT_ERROR";

        case NPP_QUAD_ERROR:
            return "NPP_QUAD_ERROR";

        case NPP_WRONG_INTERSECTION_ROI_ERROR:
            return "NPP_WRONG_INTERSECTION_ROI_ERROR";

        case NPP_NOT_EVEN_STEP_ERROR:
            return "NPP_NOT_EVEN_STEP_ERROR";

        case NPP_INTERPOLATION_ERROR:
            return "NPP_INTERPOLATION_ERROR";

        case NPP_RESIZE_FACTOR_ERROR:
            return "NPP_RESIZE_FACTOR_ERROR";

        case NPP_HAAR_CLASSIFIER_PIXEL_MATCH_ERROR:
            return "NPP_HAAR_CLASSIFIER_PIXEL_MATCH_ERROR";

        case NPP_MEMFREE_ERR:
            return "NPP_MEMFREE_ERR";

        case NPP_MEMSET_ERR:
            return "NPP_MEMSET_ERR";

        case NPP_MEMCPY_ERROR:
            return "NPP_MEMCPY_ERROR";

        case NPP_MEM_ALLOC_ERR:
            return "NPP_MEM_ALLOC_ERR";

        case NPP_HISTO_NUMBER_OF_LEVELS_ERROR:
            return "NPP_HISTO_NUMBER_OF_LEVELS_ERROR";

        case NPP_MIRROR_FLIP_ERR:
            return "NPP_MIRROR_FLIP_ERR";

        case NPP_INVALID_INPUT:
            return "NPP_INVALID_INPUT";

        case NPP_ALIGNMENT_ERROR:
            return "NPP_ALIGNMENT_ERROR";

        case NPP_STEP_ERROR:
            return "NPP_STEP_ERROR";

        case NPP_SIZE_ERROR:
            return "NPP_SIZE_ERROR";

        case NPP_POINTER_ERROR:
            return "NPP_POINTER_ERROR";

        case NPP_NULL_POINTER_ERROR:
            return "NPP_NULL_POINTER_ERROR";

        case NPP_CUDA_KERNEL_EXECUTION_ERROR:
            return "NPP_CUDA_KERNEL_EXECUTION_ERROR";

        case NPP_NOT_IMPLEMENTED_ERROR:
            return "NPP_NOT_IMPLEMENTED_ERROR";

        case NPP_ERROR:
            return "NPP_ERROR";

        case NPP_SUCCESS:
            return "NPP_SUCCESS";

        case NPP_WARNING:
            return "NPP_WARNING";

        case NPP_WRONG_INTERSECTION_QUAD_WARNING:
            return "NPP_WRONG_INTERSECTION_QUAD_WARNING";

        case NPP_MISALIGNED_DST_ROI_WARNING:
            return "NPP_MISALIGNED_DST_ROI_WARNING";

        case NPP_AFFINE_QUAD_INCORRECT_WARNING:
            return "NPP_AFFINE_QUAD_INCORRECT_WARNING";

        case NPP_DOUBLE_SIZE_WARNING:
            return "NPP_DOUBLE_SIZE_WARNING";

        case NPP_ODD_ROI_WARNING:
            return "NPP_ODD_ROI_WARNING";

        case NPP_WRONG_INTERSECTION_ROI_WARNING:
            return "NPP_WRONG_INTERSECTION_ROI_WARNING";
    }

    return "<unknown>";
}
#endif

template< typename T >
bool check(T result, char const *const func, const char *const file, int const line)
{
    if (result)
    {
        fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n",
                file, line, static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func);
        /*
                std::stringstream ss;
                std::string msg("CUDA error at ");
                msg += file;
                msg += ":";
                ss << line;
                msg += ss.str();
                msg += " code=";
                ss << static_cast<unsigned int>(result);
                msg += ss.str();
                msg += " (";
                msg += _cudaGetErrorEnum(result);
                msg += ") \"";
                msg += func;
                msg += "\"";
                //throw msg;
                std::cerr  << msg <<"\n";
        */
        return true;
    }
    else
    {
        return false;
    }
}

#ifdef __DRIVER_TYPES_H__
// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(val)           check ( (val), #val, __FILE__, __LINE__ )

// This will output the proper error string when calling hipGetLastError
#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)

inline void __getLastCudaError(const char *errorMessage, const char *file, const int line)
{
    hipError_t err = hipGetLastError();

    if (hipSuccess != err)
    {
        fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
                file, line, errorMessage, (int)err, hipGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}
#endif

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif

// Beginning of GPU Architecture definitions
inline int _ConvertSMVer2Cores(int major, int minor)
{
    // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
    typedef struct
    {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] =
    {
        { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
        { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
        { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
        { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
        {   -1, -1 }
    };

    int index = 0;

    while (nGpuArchCoresPerSM[index].SM != -1)
    {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
        {
            return nGpuArchCoresPerSM[index].Cores;
        }

        index++;
    }

    // If we don't find the values, we default use the previous one to run properly
    printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[7].Cores);
    return nGpuArchCoresPerSM[7].Cores;
}
// end of GPU Architecture definitions

#ifdef __CUDA_RUNTIME_H__
// General GPU Device CUDA Initialization
inline int gpuDeviceInit(int devID)
{
    int deviceCount;
    checkCudaErrors(hipGetDeviceCount(&deviceCount));

    if (deviceCount == 0)
    {
        fprintf(stderr, "gpuDeviceInit() CUDA error: no devices supporting CUDA.\n");
        exit(EXIT_FAILURE);
    }

    if (devID < 0)
    {
        devID = 0;
    }

    if (devID > deviceCount-1)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, ">> %d CUDA capable GPU device(s) detected. <<\n", deviceCount);
        fprintf(stderr, ">> gpuDeviceInit (-device=%d) is not a valid GPU device. <<\n", devID);
        fprintf(stderr, "\n");
        return -devID;
    }

    hipDeviceProp_t deviceProp;
    checkCudaErrors(hipGetDeviceProperties(&deviceProp, devID));

    if (deviceProp.computeMode == hipComputeModeProhibited)
    {
        fprintf(stderr, "Error: device is running in <Compute Mode Prohibited>, no threads can use ::hipSetDevice().\n");
        return -1;
    }

    if (deviceProp.major < 1)
    {
        fprintf(stderr, "gpuDeviceInit(): GPU device does not support CUDA.\n");
        exit(EXIT_FAILURE);
    }

    checkCudaErrors(hipSetDevice(devID));
    printf("gpuDeviceInit() CUDA Device [%d]: \"%s\n", devID, deviceProp.name);

    return devID;
}

// This function returns the best GPU (with maximum GFLOPS)
inline int gpuGetMaxGflopsDeviceId()
{
    int current_device     = 0, sm_per_multiproc  = 0;
    int max_compute_perf   = 0, max_perf_device   = 0;
    int device_count       = 0, best_SM_arch      = 0;
    hipDeviceProp_t deviceProp;
    hipGetDeviceCount(&device_count);

    // Find the best major SM Architecture GPU device
    while (current_device < device_count)
    {
        hipGetDeviceProperties(&deviceProp, current_device);

        // If this GPU is not running on Compute Mode prohibited, then we can add it to the list
        if (deviceProp.computeMode != hipComputeModeProhibited)
        {
            if (deviceProp.major > 0 && deviceProp.major < 9999)
            {
                best_SM_arch = MAX(best_SM_arch, deviceProp.major);
            }
        }

        current_device++;
    }

    // Find the best CUDA capable GPU device
    current_device = 0;

    while (current_device < device_count)
    {
        hipGetDeviceProperties(&deviceProp, current_device);

        // If this GPU is not running on Compute Mode prohibited, then we can add it to the list
        if (deviceProp.computeMode != hipComputeModeProhibited)
        {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
            {
                sm_per_multiproc = 1;
            }
            else
            {
                sm_per_multiproc = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor);
            }

            int compute_perf  = deviceProp.multiProcessorCount * sm_per_multiproc * deviceProp.clockRate;

            if (compute_perf  > max_compute_perf)
            {
                // If we find GPU with SM major > 2, search only these
                if (best_SM_arch > 2)
                {
                    // If our device==dest_SM_arch, choose this, or else pass
                    if (deviceProp.major == best_SM_arch)
                    {
                        max_compute_perf  = compute_perf;
                        max_perf_device   = current_device;
                    }
                }
                else
                {
                    max_compute_perf  = compute_perf;
                    max_perf_device   = current_device;
                }
            }
        }

        ++current_device;
    }

    return max_perf_device;
}


// Initialization code to find the best CUDA Device
inline int findCudaDevice(int argc, const char **argv)
{
    hipDeviceProp_t deviceProp;
    int devID = 0;

    // If the command-line has a device number specified, use it
    if (checkCmdLineFlag(argc, argv, "device"))
    {
        devID = getCmdLineArgumentInt(argc, argv, "device=");

        if (devID < 0)
        {
            printf("Invalid command line parameter\n ");
            exit(EXIT_FAILURE);
        }
        else
        {
            devID = gpuDeviceInit(devID);

            if (devID < 0)
            {
                printf("exiting...\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    else
    {
        // Otherwise pick the device with highest Gflops/s
        devID = gpuGetMaxGflopsDeviceId();
        checkCudaErrors(hipSetDevice(devID));
        checkCudaErrors(hipGetDeviceProperties(&deviceProp, devID));
        printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);
    }

    return devID;
}

// General check for CUDA GPU SM Capabilities
inline bool checkCudaCapabilities(int major_version, int minor_version)
{
    hipDeviceProp_t deviceProp;
    deviceProp.major = 0;
    deviceProp.minor = 0;
    int dev;

    checkCudaErrors(hipGetDevice(&dev));
    checkCudaErrors(hipGetDeviceProperties(&deviceProp, dev));

    if ((deviceProp.major > major_version) ||
        (deviceProp.major == major_version && deviceProp.minor >= minor_version))
    {
        printf("> Device %d: <%16s >, Compute SM %d.%d detected\n", dev, deviceProp.name, deviceProp.major, deviceProp.minor);
        return true;
    }
    else
    {
        printf("No GPU device was found that can support CUDA compute capability %d.%d.\n", major_version, minor_version);
        return false;
    }
}
#endif

// end of CUDA Helper Functions


#endif
