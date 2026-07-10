#pragma once
#include <cassert>
template <typename T> struct cudaPointer {
  T *dev_pointer = nullptr; T *host_pointer = nullptr; int size = 0;
  void allocate(int count) { size = count; checkHipErrors(hipMalloc(reinterpret_cast<void **>(&dev_pointer), size*sizeof(T))); checkHipErrors(hipHostMalloc(reinterpret_cast<void **>(&host_pointer), size*sizeof(T))); assert(dev_pointer && host_pointer); }
  void free() { if (dev_pointer) checkHipErrors(hipFree(dev_pointer)); if (host_pointer) checkHipErrors(hipHostFree(host_pointer)); dev_pointer=nullptr; host_pointer=nullptr; size=0; }
  void htod(int count) { checkHipErrors(hipMemcpy(dev_pointer, host_pointer, count*sizeof(T), hipMemcpyHostToDevice)); }
  void htod() { htod(size); }
  void dtoh(int count) { checkHipErrors(hipMemcpy(host_pointer, dev_pointer, count*sizeof(T), hipMemcpyDeviceToHost)); }
  void dtoh() { dtoh(size); }
  T &operator[](int i) { return host_pointer[i]; }
  operator T *() { return dev_pointer; }
};
