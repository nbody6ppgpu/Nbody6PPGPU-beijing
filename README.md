[![autotest status](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/actions/workflows/autotest.yml/badge.svg)](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/actions/workflows/autotest.yml)
[![Paper](https://badgen.net/badge/NASA%20ads/1999PASP..111.1333A/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/1999PASP..111.1333A/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/1999JCoAM.109..407S/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..407S/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2005MNRAS.363..293H/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.363..293H/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2012MNRAS.424..545N/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.424..545N/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2015MNRAS.450.4070W/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.4070W/abstract)
[![Paper](https://badgen.net/badge/NASA%20ads/2022MNRAS.511.4060K/blue?icon=https://ui.adsabs.harvard.edu/styles/img/transparent_logo.svg)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4060K/abstract)
<!-- [![Paper](https://badgen.net/badge/arXiv/0000.0000/green?icon=https://static.arxiv.org/static/browse/0.3.4/images/arxiv-logo-one-color-white.svg )](https://arxiv.org/abs/xxxxx) -->

This is Nbody6++GPU - Beijing version, an N-body star cluster simulation code, maintained by Rainer Spurzem (spurzem@nao.cas.cn) and team. 

The code is an offspring of Sverre Aarseth's direct N-body codes see www.sverre.com . 

This is the code suitable for parallel and GPU accelerated runs on supercomputers and workstations. Before we give some more practical help, please read the following disambiguation; there is another github of Nbody6++GPU:

LW: https://github.com/nbodyx/  - if interested please contact and collaborate
with Long Wang longwang.astro@live.com
RS: https://github.com/nbody6ppgpu - if interested please contact and
collaborate with Rainer Spurzem spurzem@ari.uni-heidelberg.de spurzem@nao.cas.cn

Here is an example of current differences between the code version (May 2023), more changes and differences may occur in the future, if in doubt, ask the authors.

1. LW: implementation of Milky Way potential following the MWPotential2014 in
Galpy (Bovy 2015).
2. RS: implementation of spin and mass dependent recoil kicks after GW merger
(Arca Sedda et al. 2023 subm. MNRAS)
3. LW:  implementation of python data reading interface for PeTar analysis tool.
4. RS: use of HDF5 output files with python data reading interfaces
5. RS: Namelist based input format, allowing also to read all stellar evolution
and binary / collision parameters.
6. LW and RS: Some bug fixes related to Roche and GR radiation, in both versions
slightly different ways.
7. LW and RS: implementation of BSE from Banerjee et al. 2019

-------------------
# Installation
## Get the code
```bash
git clone git@github.com:nbody6ppgpu/Nbody6PPGPU-beijing
```
1. This downloads the `stable` branch. The `stable` branch include major versions, and the `dev` branch include the most recent updates and bugfix. Changes in `dev` branch are merged to `stable` regularly.
2. If you want the most recent version, use 
``` bash
git clone -b dev git@github.com:nbody6ppgpu/Nbody6PPGPU-beijing
```
or run `git switch dev` after you `clone` without `-b dev` param. 

## Configure for compile

```bash
./configure [options]
```
0. TL;DR: to quickly start on your personal computer, you may use `./configure --enable-mcmodel=large --with-par=b1m --disable-gpu --disable-mpi`, and jump to the next section [Compile the code](#Compile-the-code)
1. We recommend using `--enable-mcmodel=large` to allows the program to use much resources. 
2. `--with-par=b1m` allows up to 1 million particle simulation. In case that your computer has very small memory (<4GB) and your star cluster has a small particle number, you may use smaller value (check ./configure --help for possible value for `--with-par`)
3. If you run NBODY6++GPU on your personal computer or workstation rather than computer clusters, MPI can be disabled by append `--disable-mpi` to the command above.
4. In the following cases, you may need to append `--disable-gpu` 
- The computer has no NVIDIA GPU
- The computer has NVIDIA GPU but did not install CUDA compiler (Test: type `nvcc --version` in your terminal. If you see information about NVIDIA compiler, then it is installed. If you see errors like "nvcc: command not found" then it is not installed)
- Your simulation has relatively small particle number (<50000). The code is for up to one million bodies with many initial binaries. In the case of small particle number, GPU can hardly boost the simulation and can sometimes slow it down.
5. You may set `--prefix=[install path]` to specify the location to install the executable.
6. HDF5 is an efficient storage scheme, which is useful during large-scale or long-time simulations to boost the simulation and save disk spaces. Once enabled, the basic particle data (mass, position, velocity) and stellar evolution data will be stored in `.h5part` files, which may need extra tools to read. HDF5 is recommended but not necessary. You need to install additional libraries to use HDF5. For example, in Debian based Linux `sudo apt-get install libhdf5-openmpi-dev libhdf5-dev`. Note that `--enable-hdf5` in configure command is NOT working well. You have to configure first without it, then edit build/Makefile (see example in Makefile.save.hdf5 , HDF5_DIR has to be defined by you or by system). Then make.
7. The configure script written by Long Wang has a multitude of further options, check with `./configure --help` or feel free to ask any question in [our discussion](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/discussions).

## Compile the code

```bash
make clean
make -j 
```

After `make` you can find the executable in `build/`, named `nbody6++.[configure-options]`, where the suffix depends on your configure option (MPI, GPU, HDF5, SIMD, etc), for example `nbody6++.avx.mpi.gpu`

If you have specified `--prefix=[install path]` during configure, you may want

```bash
make install
```
and add the installation path to your `$PATH` environment variable.

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

    By default, the program uses all CPU threads (which is usually 2 × number of CPU cores). For better performance, `OMP_NUM_THREADS` should not be too large, and cannot go beyond 32. In case you want to use fewer threads, especially when your computer has more than 32 cores (per node), you need to restrict `OMP_NUM_THREADS` 
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

# Documentation
To understand the diagnostic information and columns of each output file, please read the documentations at
https://www.overleaf.com/read/hcmxcyffjkzq

You are also welcomed to ask any question in [our discussion](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/discussions)

# Data analysis
Some Jupyter notebooks for simple data analysis are provided in [examples/](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/tree/stable/examples). You can check [the readme file there](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/tree/stable/examples) to get started.

# Tips
 - Before a simulation, it is always recommended to set `ulimit -s unlimited` before the simulation to avoid segmentation fault. 

 - The environment variable OMP_NUM_THREADS has to be set to the desired value of OpenMP threads per MPI process. (Maybe your system has it predefined). I also recommend to set
 OMP_STACKSIZE=4096M the shell where you run the code.

 - It is inefficient (and even more error prone) for particle numbers below about 50k-100k particles (depending on hardware). For smaller N you are advised to disable GPU, or use Nbody6 and Nbody6GPU for single node/process.
 
 - It is recommended to provide a dat.10 file in N-body input format (see manual). Such file can be produced by other programs, like McLuster.

# Seleted References:
 - https://ui.adsabs.harvard.edu/abs/1999PASP..111.1333A/abstract (Aarseth: NBODY1 to NBODY6)
 - https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..407S/abstract (Spurzem on NBODY6++)
 - https://ui.adsabs.harvard.edu/abs/2005MNRAS.363..293H/abstract (Hurley+ on SSE/BSE, earlier references therein)
 - https://ui.adsabs.harvard.edu/abs/2012MNRAS.424..545N/abstract (Nitadori+: NBODY6GPU)
 - https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.4070W/abstract (Wang+: NBODY6++GPU)
 - https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4060K/abstract (Kamlah+: More on current stellar evol.)

# For contributors
```git clone -b dev git@github.com:nbody6ppgpu/Nbody6PPGPU-beijing```
It would automatically switch to dev branch after downloading.

Sources are in `src/Main/`. Due to urgent bug fixes few routines are later than Dec2020. 

Git system does not preserve the modification time of files, but the modification time of some ancient files (created before this project was brought to Git) may be valuable information for developers. If you need this info, run `python3 restore_mtime.py` after `git clone` and each `git pull`. It will `touch` each file with their real last modification time.

# Known Problems:
 1. For systems with more than one GPU on one node the association of MPI rank id and GPU bus id is not 
      well defined, will be improved in next version.
 2. Runs with a million or more bodies and huge numbers of binaries (5% or more) use extreme amounts of
      computing time for the KS binaries (much much more than should be expected). We work on this. 
 3. On some systems heap and stack management when using OpenMP and MPI together seem to produce very
      strange errors and segmentation faults. The exact reason is not known; we work on this.
 4. Currently using standard OpenMP WITHOUT sse or avx does not work. (it means for configure --disable-simd , but --enable-omp). It uses routines nbint.F instead of special sse or avx routines for neighbour force. We are working on that.
 
 5. Using much more than one million particles (up to ten million) is still not fully supported. configure already allows --with-par=4m  , 8m, 10m, b4m, b8m, b10m . Runs of that size may still fail, depending on your hardware and software environment; also the code may still have some glitches (wrong printout, insufficient vector space allocation);  test and work is ongoing.
 
 6. Many stellar evolution and other parameters are still compiled into the code (see Table A1 in Kamlah et al. 2022, and parameter FctorCl in Rizzuto et al. 2021), mxns0,1 masses of neutron stars; it is the responsibility of the user to keep them all consistent at compile time (for example  mxns and FctorCl are defined in two routines independently, see hrplot, coal, mix). We are working to prepare a nice Fortran NAMELIST style input for ALL parameters (the ones from the current input file, and the ones currently compiled in). That will work like in the style of an .ini file with "key=value" pairs and default values.

 7. Currently the use of KZ(7) ge 4 is not working; it produces wrong SIGR2, SIGT2, VROT both in output file and lagr.7. KZ(7) le 3 is ok.
    
 8. For the latest Aug2022 stable version, there may be problems when using some extreme initial conditions (e.g. multiple massive black holes). We work on this.

# Disclaimer
 This code and the documentation is given without warranty, hopefully it is helpful. All may contain errors.
