#pragma once
#include <cassert>
#include <sys/time.h>
#define NAN_CHECK(value) assert((value) == (value))
#define PROFILE
static double get_wtime() { timeval tv; gettimeofday(&tv, nullptr); return tv.tv_sec + 1.e-6 * tv.tv_usec; }
