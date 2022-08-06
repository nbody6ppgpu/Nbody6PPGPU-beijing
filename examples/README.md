# Examples

The example directory contains:

 - 4 Jupyter notebooks which can be used as reference on how to read some nbody6++gpu output. See [below](#jupyter-notebooks) for more information.
   1. [Basics](https://github.com/kaiwu-astro/Nbody6PPGPU-beijing/blob/stable/examples/01_Basics.ipynb): Reads data from the standard output file,
   2. [Hertzsprung-Russell-diagram](https://github.com/kaiwu-astro/Nbody6PPGPU-beijing/blob/stable/examples/02_Hertzsprung%E2%80%93Russell_diagram.ipynb): Reading data from another output file, here: read data from sev.83_[...] files for an HRD.
   3. [HDF5 Basics](https://github.com/kaiwu-astro/Nbody6PPGPU-beijing/blob/stable/examples/03_HDF5_Basics.ipynb): Intro on how to use the `snap.40_[...]` files, reading multiple files to create time dependent plots (here: energy as example), reading a single file for an HRD. 
   4. [HDF5 Single file](https://github.com/kaiwu-astro/Nbody6PPGPU-beijing/blob/dev/examples/readhdf5.ipynb): Reading one snapshot as an example, and creating an N-body style enviroment with variables like ZMBAR, VSTAR, BODY, NAME, X1,... V1, ... ASPN, ... and so on. Includes three example plots (phase space, z-angular momentum, HRD).
   
   The notebooks can also be found as an html version in the `jupyter_html` subdirectory

 - example_1k: prepared folder for an example run with 1000 particles
 - input_files: Sample input files for a run with 100,000 particles including the `dat.10` file with initial conditions.



## Jupyter notebooks
### Preparations
Before examining the jupyter notebooks, the `example_1k` run should be done. In order to do the run, nbody must be build first.

#### Building
Make sure you are in the nbody git root directory and **not** in the `examples` folder.

#### Configuring hdf5
*Only needed for the hdf5 jupyter notebooks!*

As unfortunately the `--enable-hdf5` flag is broken in the configure script, you need to figure out, how to include this yourself. Both the header files path needs to be added to the `HDF5_FLAGS` using the -I switch, as well as the path to the shared library must be added via the `-L` switch in the `build/Makefile`. The `lhdf5_fortran` flag should be the same on every system. Then make sure to add the `HDF5_FLAGS` to the Fortran compiler flags `FFLAGS`.

E.g. on kepler, the `build/Makefile` should read

    HDF5_FLAGS = -D H5OUTPUT -I/usr/include/openmpi-x86_64/ -L/usr/lib64/openmpi/lib/ -lhdf5_fortran
	[...]
	FFLAGS = [...] ${HDF5_FLAGS}

#### Compile

**Note:** Both flags `--disable-mpi` `--disable-gpu` should be used for testing purposes only (e.g. the small example_1k run used as demo here)! This will slow down the simulations in general a lot! Though they can be used for the N=1000 test run without loosing too much time, it will be incredibly slow on e.g. the N100k examples! Both flags can be safely removed on any production system, and in general it is strongly advised to do so.

    ./configure --enable-simd=sse --enable-mcmodel=large --disable-mpi --disable-gpu
	# Build 
	make clean && make -j

#### prepare example_1k run
Copy the program just built to the `examples/example_1k` directory.

Note that the file name `nbody6++.sse.hdf5` might differ, as by default it depends on the flags that were set during `configure` (e.g. it might be `nbody6++.sse.gpu.mpi.hdf5`).

	# Binary name might differ!
    cp build/nbody6++.sse.hdf5 examples/example_1k

#### start the example_1k run
Move to the `examples/example_1k` directory.

    cd examples/example_1k

Start the run. Note again that the file name `nbody6++.sse.hdf5` might differ, due to the configuration flags!

	# Binary name might differ!
    ./nbody6++.sse.hdf5 < 1k.inp 1> 1k.out 2> 1k.err

With `< 1k.inp` the input configuration is given to the nbody program, wheras with `1> 1k.out` the standard output (stdout) of the program is redirected into a file called `1k.out` in this example. Similarly with `2> 1k.err` the standard error (stderr) is redirected into the file `1k.err`.

If the program runs successfully, you can now checkout the Jupyter notebooks!
