# Pulsar registry lifecycle regressions

These tests cover the evolution half of registry retirement, including the
stable-RLOF call site that can bypass the supernova branch in `roche.f`.
Run them using a clean GNU serial build (little-endian, 4-byte record markers):

```bash
./configure --with-par=1k --enable-mcmodel=large --disable-gpu \
  --disable-mpi --disable-hdf5 --enable-simd=no
make -C build -j4
python3 tests/pulsar_lifecycle/run.py --output /tmp/pulsar-lifecycle-results
```

The output directory must not exist. `--executable` selects the integrator;
`--fc` must match the compiler used for `build/*.o`. The runner links two
unit drivers and a test-only restart generator against those objects. It
never edits source, submits scheduler jobs or modifies an existing run.
Use a fresh build directory when changing configuration. The top-level
Makefile can consider an existing executable up to date without checking
changed objects; `make -C build` checks their dependencies directly.

The checks are:

- NS-to-BH retirement, repeated calls, unchanged NSs, disabled mode,
  compaction, event identity and unchanged `IDUM1`/full `RAN2` state.
- Kai's existing physics/merger invariants, including all compacted arrays
  and retirement event fields.
- A real small-N restart integration through stable RLOF and `HRDIAG`,
  ending with a black hole, zero NS slots, and exactly one code-6 event.
- A further restart retaining that result without recreating pulsar state.
- Missing/truncated/invalid/duplicate/stale restart records, missing live-NS
  history, off-to-on restarts, and invalid `KZ(19)` combinations.
- Fresh type-13 input rejection with pulsars enabled; the same evolved
  stellar input remains accepted with pulsars disabled.
- The four portable CE/collision fixtures, checked against physical stars
  in every checkpoint, with birth logs and collision/CE event paths.
  `provisional_ns_cleanup` is reported only as an environment-sensitive
  smoke test. This suite does not claim EXPEL2 retry coverage.

## Synthetic accretion fixture

`rlof_ns_to_bh` contains 50 stars: a nearly circular binary with
`a=2.7 Rsun`, `e=0.001`, a `(2.5 - 1e-10) Msun` type-13 accretor and a
`1 Msun` main-sequence donor, plus 48 `0.3 Msun` spectators. The stellar
reference masses/ages are `10 Msun`/`1 Myr` for the NS and
`1 Msun`/`1000 Myr` for the donor. `RBAR=0.01 pc`; natal kicks and CE
pulsar accretion are disabled. The accretor starts just below the default
`MXNS=2.5 Msun` so the first mass-transfer episode crosses the threshold.
This is a regression fixture, not a model for population inference.

`seed_restart.f` runs the ordinary initialization with pulsars disabled,
then explicitly supplies test pulsar state (`P=0.05 s`, `B=1e8 G`,
`Pdot=1e-15`, registry ages zero), enables pulsars and writes a normal dump.
The helper suppresses birth logging during this synthetic initialization;
no physical formation history is claimed. This helper is never linked
into the production executable and does not add support for fresh NS ICs.
The runner copies the saved file to `comm.1` and uses the unmodified
production entry point to integrate for `0.2` N-body time units. Dump and
integration logs remain under the requested output directory, together
with `report.json`, the executable hash and compiler identity.

The restart audit deliberately retains unmatched names (escapers or hidden
subsystems) and warns. Unit tests exercise that policy. It rejects definite
inconsistencies without inventing or deleting unresolved history.

## MPI checks

Build MPI separately with the same particle limits and compiler. The serial
suite's synthetic dumps can then be used to verify one-rank MPI integration
and two-rank rejection of invalid restarts:

```bash
python3 tests/pulsar_lifecycle/check_mpi.py \
  --executable /path/to/mpi/build/nbody6++.mpi \
  --serial-results /tmp/pulsar-lifecycle-results \
  --output /tmp/pulsar-lifecycle-mpi-results
```

`--launcher` accepts additional MPI launcher arguments. These tests have
been checked with GNU Fortran 11.3.0 and OpenMPI 4.1.4. They do not claim
multi-rank integration coverage: restarting this serial-generated fixture
with two ranks and the no-SIMD MPI build hit `MPI_ERR_TRUNCATE` in both the
unmodified owner baseline `8490a2c` and the patched code. The positive
integration checks therefore use one MPI rank; all eight negative restart
checks use two ranks and must terminate with the intended diagnostic,
without hanging.
