This is Nbody6++GPU - Beijing version, an N-body star cluster simulation code, maintained by Rainer Spurzem (spurzem/at/nao.cas.cn) and team. 

The code is an offspring of Sverre Aarseth's direct N-body codes see www.sverre.com . 

This is the code suitable for parallel and GPU accelerated runs on supercomputers and workstations. 

# For contributors
```git clone -b dev git@github.com:nbody6ppgpu/Nbody6PPGPU-beijing```
It would automatically switch to dev branch after downloading.

Sources are in `src/Main/`. Due to urgent bug fixes few routines are later than Dec2020. 

Git system does not preserve the modification time of files, but the modification time of some ancient files (created before this project was brought to Git) may be valuable information for developers. If you need this info, run `python3 restore_mtime.py` after `git clone` and each `git pull`. It will `touch` each file with their real last modification time.

# For users
## Installation
```bash
git clone git@github.com:nbody6ppgpu/Nbody6PPGPU-beijing
```
This downloads the `stable` branch. The `stable` branch include major versions, and the `dev` branch include the most recent updates and bugfix. Changes in `dev` branch are merged to `stable` regularly.

If you want the most recent version (may contain bugs), use `git clone -b dev git@github.com:nbody6ppgpu/Nbody6PPGPU-beijing`, or run `git switch dev` after you `clone` without `-b dev` param. 

```bash
 ./configure --with-par=b1m --enable-simd=sse --enable-mcmodel=large 
 make clean ; make -j 
```

After make you find the executable and object files in `build/`, named `nbody6++.sse`. The suffix may change with different compilation options.

Copy the executable to the simulation directory you want

```bash
cp `ls build/nbody6++*` [your_simulation_path]
```

It is for up to one million bodies with many initial binaries. The configure script written by Long Wang has a multitude of further options, check with `./configure --help` or feel free to ask any question in [our discussion](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/discussions).

## Additional installation options

HDF5 is recommended but not necessary. To use HDF5, make sure it is installed in your computer. For example, in Debian based Linux,
```bash
apt-get install libhdf5-openmpi-dev
apt-get install libhdf5-dev
```

after that, append `--enable-hdf5` in configure command

# Documentation
https://www.overleaf.com/read/hcmxcyffjkzq

# Tips
 - The environment variable OMP_NUM_THREADS has to be set to the desired value of OpenMP threads per MPI process. (Maybe your system has it predefined). I also recommend to set
 OMP_STACKSIZE=4096M the shell where you run the code.

 - It is inefficient (and even more error prone) for particle numbers below about 50k-100k particles (depending on hardware). For smaller N you are advised to use Nbody6 and Nbody6GPU for single node/process.
 
 - It is recommended to provide a dat.10 file in N-body input format (see manual). Such file can be produced by other programs, like mcluster.


# Seleted References:
 - https://ui.adsabs.harvard.edu/abs/1999PASP..111.1333A/abstract (Aarseth: NBODY1 to NBODY6)
 - https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..407S/abstract (Spurzem on NBODY6++)
 - https://ui.adsabs.harvard.edu/abs/2005MNRAS.363..293H/abstract (Hurley+ on SSE/BSE, earlier references therein)
 - https://ui.adsabs.harvard.edu/abs/2012MNRAS.424..545N/abstract (Nitadori+: NBODY6GPU)
 - https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.4070W/abstract (Wang+: NBODY6++GPU)
 - https://ui.adsabs.harvard.edu/abs/2021arXiv210508067K/abstract (Kamlah+: More on current stellar evol.)


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
