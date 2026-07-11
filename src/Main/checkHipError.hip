#pragma once
#include <cstdio>
#include <cstdlib>
#define checkHipErrors(call) do { const hipError_t e__ = (call); if (e__ != hipSuccess) { std::fprintf(stderr, "HIP error at %s:%d: %s\n", __FILE__, __LINE__, hipGetErrorString(e__)); std::exit(EXIT_FAILURE); } } while (0)
#define checkHipLaunch() checkHipErrors(hipGetLastError())
