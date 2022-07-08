This is Nbody6++GPU - Beijing version, an N-body star cluster simulation code, maintained by Rainer Spurzem (spurzem/at/nao.cas.cn) and team. The starting commits of this repo put all old versions of Nbody6++GPU from Dec 2017 into Git system for version control.

The code is an offspring of Sverre Aarseth's direct N-body codes see www.sverre.com . 

This is the code suitable for parallel and GPU accelerated runs on supercomputers and workstations. 

# About this repo
**The `stable` branch include major versions, and the `dev` branch include the most recent updates and bugfix. Changes in `dev` branch are merged to `stable` regularly.**

**If you want the most recent version, run `git switch dev` after you `git clone` the code**

# Documentation
Manual in `doc/nbody6++_manual.pdf`


# Recommended usage:
## (optional) to enable HDF5
```bash
apt-get install libhdf5-openmpi-dev
apt-get install libhdf5-dev
```
## build
```bash
 ./configure --with-par=b1m --enable-simd=sse --enable-mcmodel=large
 make clean ; make -j 
```
 It is for up to one million bodies with many initial binaries. The configure script written by 
 Long Wang has a multitude of further options, check in it or ask.

 Sources are in src/Main/ . Due to urgent bugfixes few routines are later then Dec2020. 

 After make you find the executable and object files in build/  .
 
# Tips
 The environment variable OMP_NUM_THREADS has to be set to the desired value of 
 OpenMP threads per MPI process. (Maybe your system has it predefined). I also recommend to set
 OMP_STACKSIZE=4096M the shell where you run the code.

 It is inefficient (and even more error prone) for particle numbers below about 50k-100k particles (depending on hardware). For smaller N you are advised to use Nbody6 and Nbody6GPU for single node/process.
 
 It is recommended to provide a dat.10 file in N-body input format (see manual). Such file can be produced by other prograns, like mcluster.


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


# Disclaimer
 This code and the documentation is given without warranty, hopefully it is helpful. All may contain errors.
