[![autotest status](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/actions/workflows/autotest.yml/badge.svg)](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/actions/workflows/autotest.yml)
[![Paper](https://badgen.net/badge/NASA%20ads/1999PASP..111.1333A/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/1999PASP..111.1333A/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/1999JCoAM.109..407S/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..407S/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2005MNRAS.363..293H/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.363..293H/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2012MNRAS.424..545N/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.424..545N/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2015MNRAS.450.4070W/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.4070W/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2022MNRAS.511.4060K/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4060K/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2023LRCA....9....3S/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2023LRCA....9....3S/abstract)


<!-- [![Paper](https://badgen.net/badge/arXiv/0000.0000/green?icon=https://static.arxiv.org/static/browse/0.3.4/images/arxiv-logo-one-color-white.svg )](https://arxiv.org/abs/xxxxx) -->

This is Nbody6++GPU - Beijing version, an N-body star cluster simulation code, maintained by Rainer Spurzem (spurzem@nao.cas.cn) and team, main developers Kai Wu (kaiwu.astro@gmail.com) and Francesco Flammni Dotti (ff2415@nyu.edu).

The code is an offspring of [Sverre Aarseth's direct N-body codes](https://people.ast.cam.ac.uk/~sverre/web/pages/nbody.htm).

Any important messages from our side to the community (e.g. updates which change downward compatibility or other major issues) see in our discussion section: https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/discussions

This is the code suitable for parallel and GPU accelerated runs on supercomputers and workstations. Before we give some more practical help, please read the following disambiguation; there is another github of Nbody6++GPU:

LW: https://github.com/nbodyx/  - if interested please contact and collaborate
with Long Wang longwang.astro@live.com
RS: https://github.com/nbody6ppgpu - if interested please contact and
collaborate with Rainer Spurzem spurzem@nao.cas.cn

Here is an example of current differences between the code version (May 2023), more changes and differences may occur in the future, if in doubt, ask the authors.

1. LW: implementation of python data reading interface for PeTar analysis tool.
2. RS: implementation of spin and mass dependent recoil kicks after GW merger (Arca Sedda et al. 2023 subm. MNRAS)
3. RS: use of HDF5 output files with python data reading interfaces
4. RS: namelist based input format, allowing also to read all stellar evolution and binary / collision parameters.
5. LW and RS: implementation of Milky Way potential following the MWPotential2014 in Galpy (Bovy 2015).
6. LW and RS: Some bug fixes related to Roche and GR radiation, in both versions slightly different ways.
7. LW and RS: implementation of BSE from Banerjee et al. 2019

-------------------
# Installation
## Get the code
```bash
git clone https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing
```
1. This downloads the `stable` branch. The `stable` branch include major versions, and the `dev` branch include the most recent updates and bugfix. Changes in `dev` branch are merged to `stable` regularly.
2. If you want the most recent version, use
``` bash
git clone -b dev https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing
```
or run `git switch dev` after you `clone` without `-b dev` param.

## Configure for compile

```bash
./configure [options]
```
0. TL;DR: choose the command for your target architecture, then jump to [Compile the code](#Compile-the-code):
   - x86_64: `./configure --with-par=b1m --enable-mcmodel=large --disable-gpu`
   - ARM64/aarch64: clone SIMDe as shown below, then use `./configure --with-par=b1m --with-simde="../simde" --enable-simd=sse --disable-gpu`
1. On x86_64, we recommend `--enable-mcmodel=large` to allow the program to use more resources. ARM64 does not support this option and defaults to `mcmodel=no`.
2. `--with-par=b1m` allows up to 1 million particle simulation. In case that your computer has very small memory (<4GB) and your star cluster has a small particle number, you may use smaller value (check ./configure --help for possible value for `--with-par`)
3. MPI should always be used during compilation. In case your computer does not have it, you can install with `sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev` in Debian based Linux. The option `--disable-mpi` should only be used for debug purpose, and not for any production run.
4. In the following cases, you may need to append `--disable-gpu`
- Your simulation has relatively small particle number (<50000). The code is for up to one million bodies with many initial binaries. In the case of small particle number, GPU can hardly boost the simulation and can sometimes slow it down.
- The computer has no supported NVIDIA CUDA or AMD ROCm GPU/toolchain.

5. GPU builds default to `--with-gpu-backend=auto`, which selects `nvcc` before `hipcc`. Use `--with-gpu-backend=cuda` or `--with-gpu-backend=hip` for a strict selection; an explicit backend never falls back. CUDA accepts `--with-cuda=PREFIX` and `NVCC=/path/to/nvcc`; standard AMD ROCm accepts `--with-hip=PREFIX` and `HIPCC=/path/to/hipcc`. If neither compiler and runtime can compile/link, configure stops and suggests `--disable-gpu`. Do not combine `--disable-gpu` with an explicit CUDA/HIP backend.
6. You may set `--prefix=[install path]` to specify the location to install the executable.
7. HDF5 is an efficient storage scheme, which is useful during large-scale or long-time simulations to boost the simulation and save disk spaces. The basic particle data (mass, position, velocity) and stellar evolution data will be stored in `.h5part` files, which may need extra tools to read. HDF5 support is auto-detected and **enabled by default**, since all current output formats are built on it and the legacy ASCII output is no longer maintained. Install the serial HDF5 Fortran development package (for example, `libhdf5-dev` on Debian-based Linux) before running `./configure`. The configure script discovers include and library paths through `h5fc`; set `H5FC=/path/to/h5fc` when the serial wrapper is not first in `PATH`. Parallel HDF5 is not supported and produces warnings during both configure and make. If `h5fc` cannot be found, configure stops with instructions to install HDF5 or, if you really do not want HDF5 output, to pass `--disable-hdf5` (not recommended).
8. The configure script written by Long Wang has a multitude of further options, check with `./configure --help` or feel free to ask any question in [our discussion](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/discussions).

### ARM64 and SIMDe

ARM64 builds use [SIMDe](https://github.com/simd-everywhere/simde) to run the SSE force kernels without changing the simulation format. SIMDe v0.8.2 is the version exercised by CI. Clone it next to this repository (not inside it):

```bash
git clone --branch v0.8.2 --depth 1 https://github.com/simd-everywhere/simde.git ../simde
./configure --with-par=b1m --with-simde="../simde" --enable-simd=sse --disable-gpu
make clean && make -j
```

If `--with-simde` is omitted on ARM64, configure automatically looks for the absolute equivalent of `../simde`. Relative and absolute paths are accepted; bare `--with-simde` uses compiler system include directories. ARM64 defaults to SSE through SIMDe, so `--enable-simd=avx`, `--without-simde` (without also passing `--enable-simd=no`), and `--enable-mcmodel=small|medium|large` fail with a corrective example instead of being downgraded. `--enable-simd=no` is accepted on both ARM64 and x86_64, but see the note below about what it disables. Only x86_64 and 64-bit ARM targets are currently supported.

`--enable-simd=no` also switches OpenMP off automatically (on both architectures), because without the SSE/AVX kernels the code falls back to the plain Fortran `nbint.F` force routine, which is not thread-safe under OpenMP (Known Problems #3 below). Configure prints a warning when this happens; the resulting build is correct but single-threaded, so use `--enable-simd=no` for debugging only, never for production runs.

## Compile the code

```bash
make clean; make -j
```

After `make` you can find the executable in `build/`, named `nbody6++.[configure-options]`. CUDA keeps the `.gpu` suffix (for example `nbody6++.avx.mpi.gpu`); HIP uses `.hip`.

If you have specified `--prefix=[install path]` during configure, you may want `make install`, and add the installation path to your `$PATH` environment variable.

# Ready for your simulation

1. (If you have done `make install` you can skip this step) Copy the executable to the simulation directory you want

    ```bash
    cp `ls build/nbody6++*` [your_simulation_dir]
    ```

2. Prepare an initial condition file. For a test run, you can find example initial conditions in `examples/input_files`.

    ```bash
    cp examples/input_files/N10k_noDat10.inp [your_simulation_dir]
    ```

    This input file let NBODY6++GPU generate a star cluster with 10000 stars with Plummer model, and simulate for only 2 Myr. You can also find `N100k.inp` and its pre-generated initial particle data `dat.10` in `examples/input_files` for a 100,000 stars, 1 Gyr simulation.

    > 💡 Starting from the stable version [May2023](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/releases/tag/May2023), NBODY6++GPU changes to a fundamentally new and more flexible method of reading input data (control data, not particle data). It uses Fortran NAMELIST input, which has a key=value format. All input data can be given in any order. If you are using a old-format input file, you can use the bash script which transform the old input file into the new one ([examples/input_files/@input-transform](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/blob/stable/examples/input_files/%40input-transform)) to transform it to the new NAMELIST format. See usage inside the script.

3. CPU and memory

    In simulations with large particle number, segmentation fault may happen. To avoid this, we recommend setting a large `OMP_STACKSIZE` and disable the memory limitation.
    ```bash
    export OMP_STACKSIZE=4096M
    ulimit -s unlimited
    ```

    By default, the program uses all CPU threads (which is usually 2 × number of CPU cores). The implementation supports up to 1024 OpenMP threads, but values above 32 are rarely efficient. Benchmark the physical cores on your node and avoid using extra hyperthreads unless measurements show a benefit.
    ```bash
    export OMP_NUM_THREADS=[N_threads]
    ```

    After running them, you may want to add the these 3 commands to your shell initial file like `~/.bashrc`.

4. Finally, run it

    ```bash
    cd [your_simulation_path]
    ```

    If you have done `make install` and add the installation path to `$PATH`, run
    ```bash
    nbody6++ < N10k_noDat10.inp
    ```

    otherwise you may have copied the executable to the simulation path, run
    ```bash
    ./[your executable filename] < N10k_noDat10.inp
    ```

# Benchmark tool

`tool/benchmark.py` prepares local or Slurm benchmark matrices and collects timing results. It requires Python 3.8 or newer; YAML configuration and CSV collection additionally use PyYAML and pandas as listed in the script metadata.

```bash
python3 tool/benchmark.py --help
python3 tool/benchmark.py --generate-example-param
python3 tool/benchmark.py -N 50k,100k --mpi-per-node=1,2 --expert
```

The default input is `examples/N1m_benchmark.inp`. For Slurm, copy `examples/example.sbatch`, add site-specific account, partition, modules, and GPU options, then pass it with `--sbatch-base-path`. Local runs are successful only when the process exits with status zero and its output contains `END RUN`.

# Documentation
To understand the diagnostic information and columns of each output file, please read the documentations at
https://nbody6ppgpu.github.io/nb6-manual-pdf/latest.pdf

Which mirrors the results from the following overleaf 
https://www.overleaf.com/read/hcmxcyffjkzq#89d2bb

Please report any issues for the manual by opening an issue https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/issues . 

You are also welcomed to ask any question in [our discussion](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/discussions)

# Data analysis
Some Jupyter notebooks for simple data analysis are provided in [examples/](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/tree/stable/examples). You can check [the readme file there](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/tree/stable/examples) to get started.

# Tips
 - Before a simulation, it is always recommended to set `ulimit -s unlimited` before the simulation to avoid segmentation fault.

 - The environment variable OMP_NUM_THREADS has to be set to the desired value of OpenMP threads per MPI process. (Maybe your system has it predefined.) We also recommend to set OMP_STACKSIZE=4096M the shell where you run the code.

 - It is inefficient (and even more error prone) for particle numbers below about 50k-100k particles (depending on hardware). For smaller N you are advised to disable GPU, or use Nbody6 and Nbody6GPU for single node/process.

 - It is recommended to provide a `dat.10` file in N-body input format (see manual). Such file can be produced by other programs, like [McLuster](https://github.com/agostinolev/mcluster).

# Seleted References:
 - https://ui.adsabs.harvard.edu/abs/1999PASP..111.1333A/abstract (Aarseth: NBODY1 to NBODY6)
 - https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..407S/abstract (Spurzem on NBODY6++)
 - https://ui.adsabs.harvard.edu/abs/2005MNRAS.363..293H/abstract (Hurley+ on SSE/BSE, earlier references therein)
 - https://ui.adsabs.harvard.edu/abs/2012MNRAS.424..545N/abstract (Nitadori+: NBODY6GPU)
 - https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.4070W/abstract (Wang+: NBODY6++GPU)
 - https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4060K/abstract (Kamlah+: More on current stellar evolution)

# For contributors
```git clone -b dev git@github.com:nbody6ppgpu/Nbody6PPGPU-beijing```

Sources are in `src/Main/`.

`include/params.h`, the top-level and build `Makefile` files, `config.log`, and `config.status` are generated by `./configure` and are intentionally ignored. Change particle-size mappings in `configure.ac`, regenerate `configure` with Autoconf 2.71, and never commit local configuration output.

Git system does not preserve the modification time of files, but the modification time of some ancient files (created before this project was brought to Git) may be valuable information for developers. If you need this info, run `python3 restore_mtime.py` after `git clone` and each `git pull`. It will `touch` each file with their real last modification time.

# Known Problems:
 0. CUDA and standard AMD ROCm compile/link configuration is covered without GPU hardware, but hardware validation is still pending for both backends: N10k reaching `END RUN` without NaN/runtime errors, single/multi-GPU execution, and `GPU_LIST`. These checks must be completed on real CUDA and ROCm systems before claiming production validation. DCU/DTK is not currently supported.

 1. For systems with more than one GPU on one node the association of MPI rank id and GPU bus id is not
      well defined, will be improved in next version.

 2. Runs with a million or more bodies and huge numbers of binaries (5% or more) use extreme amounts of
      computing time for the KS binaries (much much more than should be expected). We work on this.

 3. Currently using standard OpenMP WITHOUT sse or avx does not work. (it means for configure --enable-simd=no, but with OpenMP). It uses routines nbint.F instead of special sse or avx routines for neighbour force. Since `configure` has no `--disable-omp` option, `--enable-simd=no` now automatically disables OpenMP as well (see [ARM64 and SIMDe](#arm64-and-simde) above), so this broken combination can no longer be built through configure. The underlying thread-safety issue in nbint.F itself is still unresolved. We are working on that.

 4. Many stellar evolution and other parameters are still compiled into the code (see Table A1 in Kamlah et al. 2022, and parameter FctorCl in Rizzuto et al. 2021), mxns0,1 masses of neutron stars; it is the responsibility of the user to keep them all consistent at compile time (for example  mxns and FctorCl are defined in two routines independently, see hrplot, coal, mix). We are working to prepare a nice Fortran NAMELIST style input for ALL parameters (the ones from the current input file, and the ones currently compiled in). That will work like in the style of an .ini file with "key=value" pairs and default values.

 <!-- 5. Using much more than one million particles (up to ten million) is still not fully supported. configure already allows --with-par=4m  , 8m, 10m, b4m, b8m, b10m . Runs of that size may still fail, depending on your hardware and software environment; also the code may still have some glitches (wrong printout, insufficient vector space allocation);  test and work is ongoing. -->
 <!-- 3. On some systems heap and stack management when using OpenMP and MPI together seem to produce very
      strange errors and segmentation faults. The exact reason is not known; we work on this. -->

# Disclaimer
 This code and the documentation is given without warranty, hopefully it is helpful. All may contain errors.
