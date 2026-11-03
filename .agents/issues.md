# Known Issues

| Issue | Status | Solution |
|------|------|----------|
| GPU assignment on multi-GPU nodes | In progress | Bind GPU bus IDs manually |
| Poor performance with binary fraction >20% | Under investigation | Reduce the number of binaries or use larger hardware |
| `--disable-simd` + `--enable-omp` does not work | In progress | Use SSE/AVX |
| Many parameters are hard-coded at compile time | Being improved | Use full Namelist input in the future |
| `KZ(7) >= 4` gives wrong output | Known | Use `KZ(7) <= 3` |
| HDF5 configure options do not take effect | Known | Edit `build/Makefile` manually |
| CUDA/HIP verification on real hardware | To be verified | Verify separately: N10k `END RUN`, no NaN/runtime error, single/multi-GPU and `GPU_LIST` |
