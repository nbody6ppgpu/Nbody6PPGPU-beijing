# Build and Configuration

## Configuration

```bash
./configure [options]
```

| Option | Description | Recommended |
|------|------|--------|
| `--enable-mcmodel=ARG` | Memory model (small/medium/large) | `large` |
| `--with-par=ARG` | Maximum particle number (1k/10k/100k/1m/b1m/b4m/b10m) | `b1m` (10 million) |
| `--disable-gpu` | Disable GPU | Only if N<50k or no GPU |
| `--with-gpu-backend=ARG` | Backend (`auto/cuda/hip`, auto prefers CUDA) | `auto` |
| `--with-cuda=PREFIX` / `--with-hip=PREFIX` | CUDA / standard ROCm root directory | Per install location |
| `--enable-simd=ARG` | SIMD optimisation (sse/avx/no); `no` automatically disables OpenMP (the nbint.F fallback path is not safe under OpenMP, see `.agents/issues.md`) | `avx` (if supported) |
| `--with-simde[=PATH]` | Use SIMDe headers (required by default on ARM64, unless `--enable-simd=no` is also passed) | Use `v0.8.2` on ARM64 |
| `--disable-hdf5` | Disable HDF5 output (auto-detected and enabled by default) | **Not recommended** (the new output formats are all HDF5-based) |
| `--disable-mpi` | Disable MPI | **Not recommended** (debugging only) |
| `--enable-debug` | Debug mode | Use during development |
| `--prefix=PATH` | Install path | Custom |

```bash
# Personal laptop test (no GPU, no MPI)
./configure --enable-mcmodel=large --with-par=b1m --disable-gpu --disable-mpi

# Typical production environment
./configure --enable-mcmodel=large --with-par=b1m

# ARM64 (the SIMDe path can also be an absolute path or a system install path)
git clone --branch v0.8.2 --depth 1 https://github.com/simd-everywhere/simde.git ../simde
./configure --with-par=b1m --with-simde=../simde --enable-simd=sse --disable-gpu
```

On ARM64 the unsupported `-mcmodel` is turned off by default, and the SSE kernel is selected by default when `--with-simde` is given; requesting AVX on ARM64 explicitly falls back to SSE. x86_64 keeps the native SSE/AVX paths.

## Regenerating configure

After changing `configure.ac` or the `macro/*.m4` files it `m4_include`s, you must run `tool/regen-configure.sh` to regenerate `configure` and commit it together.
`configure` is committed to the repo (cluster users without autoconf installed must still be able to run `./configure`); CI regenerates it and requires a byte-for-byte match with the committed version (`tool/regen-configure.sh --check`).
Do not run `autoconf` directly: the output of autoconf 2.71 differs between distributions (e.g. `ac_cv_prog_cxx_11` vs `ac_cv_prog_cxx_cxx11`).
The script runs in a container; the base image is pinned by digest, autoconf and m4 are pinned by version and checked at run time, and CI uses the same script.
It is not fully hermetic (the apt index is still mutable), but once the toolchain drifts it fails with an error instead of producing an unreadable diff.
Also note that `dnl` swallows only the text from itself to the end of the line, not the whitespace before it, so `dnl` comments in `configure.ac` must start at column 1.

## Compile and install

`make clean && make -j`; `make install` to install.

Executables are in `build/`, named `nbody6++.[simd].[mpi].[gpu].[hdf5]`, e.g. `nbody6++.avx.mpi.gpu`.

CUDA builds get the suffix `.gpu`, HIP builds `.hip`. `NVCC` and `HIPCC` can override the compilers; an explicit backend does not fall back, and `--disable-gpu` must not be combined with an explicit `cuda/hip`. HIP currently targets standard AMD ROCm only; N10k, single/multi-GPU and `GPU_LIST` on real CUDA/ROCm hardware are yet to be verified.

## CI

`.github/workflows/`: `autotest.yml` (automatic test build on push/PR, including the byte-for-byte `configure` check), `arm_test.yml` (simulation test on a native ARM64 runner), `release_on_tag.yml` (release on tag push), `sync-stardisk.yml` (scheduled/manual sync of specific branches).
