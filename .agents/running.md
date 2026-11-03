# Running, Output, Debugging and Tuning

## Run example

```bash
export OMP_STACKSIZE=4096M
export OMP_NUM_THREADS=16  # Implementation limit is 1024; usually keep it at or below 32 and tune by measurement; avoid using all logical cores (hyper-threading)
ulimit -s unlimited

cp examples/input_files/N10k_noDat10.inp ./
./build/nbody6++.avx.mpi.gpu < N10k_noDat10.inp
```

On out-of-memory: `--enable-mcmodel=large` + the `OMP_STACKSIZE` above + `ulimit -s unlimited`.

## Example inputs (`examples/input_files/`)

| File | Particles | Simulated time | Purpose |
|------|--------|----------|------|
| `N10k_noDat10.inp` | 10,000 | 2 Myr | Quick test (Plummer model) |
| `N100k.inp` | 100,000 | 1 Gyr | Production run example |
| `dat.10` | - | - | Pre-generated initial conditions |

## Output files

| File | Content |
|------|------|
| `OUT3` | Main output file (orbits, energy, collisions, etc.) |
| `snap.40_*` | HDF5 snapshot files |
| `sev.83_*` | Stellar evolution data |
| `lagr.7` | Lagrangian radii |
| `bev.82` | Binary evolution events |

## Data analysis

Jupyter notebooks under `examples/`: `01_Basics.ipynb` (basic reading), `02_Hertzsprung–Russell_diagram.ipynb` (HR diagram), `03_HDF5_Basics.ipynb` (HDF5 handling), `readhdf5.ipynb` (reading a single HDF5 file).

## Debugging

```bash
# Fortran
./configure --enable-debug --disable-gpu --disable-mpi && make clean && make -j
gdb build/nbody6++.debug
(gdb) run < examples/input_files/N10k_noDat10.inp

# GPU
cuda-memcheck ./build/nbody6++.avx.mpi.gpu < input.inp
cuda-gdb ./build/nbody6++.avx.mpi.gpu
```

## Performance tuning

- Vectorise with SSE/AVX; put the regular force computation on the GPU; use `!$OMP PARALLEL DO` in force loops
- Optimise the neighbour search radius; tune `DTMIN`, `DTMAX`
- Disable the GPU when N < 50k (`--disable-gpu`); multi-GPU nodes need manual binding
- Performance drops significantly when the binary fraction is >5%; consider tuning `KMAX`
