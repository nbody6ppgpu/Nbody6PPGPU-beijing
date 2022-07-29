# Examples



## Jupyter notebooks
In this folder are some jupyter notebooks on how to read the output data, and generate some intereseting plots.
### Preparations
Before examining the jupyter notebooks, the `example_1k` run should be done. In order to do the run nbody must be build first.

#### Building
Make sure you are in the nbody git root directory and **not** in the `examples` folder.

For the snap part of the jupyter notebooks, the `--enable-hdf5` flag is important! The other configure flags may be exchanged/removed.

    ./configure --with-par=b1m --enable-simd=sse --enable-mcmodel=large --disable-mpi --disable-gpu --enable-hdf5
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

With `< 1k.inp` the input configuration is given to the nbody program, wheras with `1> 1k.out` the standard output (stdout) of the program is redirected into a file in this example called `1k.out`. Similarly with `2> 1k.err` the standard error (stderr) is redirected into the file `1k.err`

If the program runs successfully, you can now checkout the Jupyter notebooks!
