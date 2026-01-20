#include <assert.h>
// #include <cutil.h>

template <typename T>
struct cudaPointer{
  T *dev_pointer;
  T *host_pointer;
  int size;
  cudaPointer(){
    dev_pointer = NULL;
    host_pointer = NULL;
    size = 0;
  }
  //  ~cudaPointer(){
  // free();
  //  }
  void allocate(int _size){
    size = _size;
    void *p;
    checkHipErrors(hipMalloc(&p, size * sizeof(T)));
    assert(p);
    dev_pointer = (T*)p;
    checkHipErrors(hipHostMalloc(&p, size * sizeof(T)));
    assert(p);
    host_pointer = (T*)p;
  }
  void free(){
    checkHipErrors(hipFree(dev_pointer));
    checkHipErrors(hipHostFree(host_pointer));
    dev_pointer = NULL;
    host_pointer = NULL;
    size = 0;
  }
  void htod(int count){
    checkHipErrors(hipMemcpy(dev_pointer, host_pointer, count * sizeof(T), hipMemcpyHostToDevice));
  }
  void htod(){
    this->htod(size);
  }
  void dtoh(int count){
    checkHipErrors(hipMemcpy(host_pointer, dev_pointer, count * sizeof(T), hipMemcpyDeviceToHost));
  }
  void dtoh(){
    this->dtoh(size);
  }
  T &operator [] (int i){
    return host_pointer[i];
  }
  operator T* (){
    return dev_pointer;
  }
};
