// https://www.cnblogs.com/zhcnfyy/p/15184193.html

#define checkHipErrors( a ) do { \
    if (hipSuccess != (a)) { \
    fprintf(stderr, "Hip runtime error in line %d of file %s \
    : %s \n", __LINE__, __FILE__, hipGetErrorString(hipGetLastError()) ); \
    /*exit(EXIT_FAILURE);*/ \
    } \
    } while(0);

