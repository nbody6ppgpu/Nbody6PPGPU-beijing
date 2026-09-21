# Pulsar registry end-to-end fixtures

This directory contains five deterministic 50-body initial conditions for the
pulsar-registry paths exercised by the common-envelope and direct-collision
changes. Each case has two or three target stars and distant 0.3-solar-mass
spectators. The repository stores inputs and expected outcomes only; generated
logs, checkpoints, executables, and Slurm files are deliberately excluded.

The fixtures were run and audited at source commit
`c1549ca2f7d92995e54cf007eaed80d37df9289a`. The validated serial executable
had SHA256
`47261c1118b071a84979806025beef17d9b0c916e91e8643c7bde145584b0015`.
`manifest.json` records the expected results and the checksums of the exact
validated inputs.

## Cases

| Case | Initial target | Validated result |
| --- | --- | --- |
| `ce_survive` | K=8 helium HG (2.368 Msun) + K=1 MS (2.909 Msun), a=18 Rsun, e=0.016082 | Two COMENV passes; a 1.277583525-Msun K=13 star and its MS companion survive. The registry contains exactly the live NS under NAME=1. |
| `ce_merge` | Same pair, a=12 Rsun | Two COMENV passes followed by COAL; one 1.277583525-Msun K=13 survivor remains. The disappeared member has no registry entry. |
| `wd_collision` | K=12 ONe WD (1.2 Msun) + K=11 CO WD (0.5 Msun) | CMBODY/MIX forms one 1.26-Msun K=13 star under NAME=1, with no duplicate registry entry. |
| `existing_ns_collision_migration` | The WD pair above, followed by a head-on 1.3-Msun K=1 intruder | The first collision forms a 1.26-Msun NS under NAME=1. The second makes a 1.91-Msun NS and moves the same registry slot to the surviving NAME=3 while preserving its spin, field, and formation age. All 101 checkpoints pass. |
| `provisional_ns_cleanup` | K=9 stripped helium giant (1.50639 Msun; reference mass 2.4064888718216153 Msun) + K=1 MS (2.909 Msun), a=8 Rsun, e=0.7045462436529022 | CHAOS enters EXPEL. The first COMENV pass produces provisional K=13; enforced CE then merges to final K=1. The temporary NS is never registered. All 81 checkpoints have NSCOUNT=0; max abs(ERRTOT)=2.23879e-9. |

For every case, `fixture.inp`, `dat.10`, and `datsev.21` must be copied
together. `KZ(22)=2` reads the supplied normalized coordinates and
`KZ(19)=4` reads the supplied evolved-star state. The runs use fixed seed
`NRAND=43532`, zero NS/WD kicks, `KZ(29)=2`, `KZ(50)=1`, and
`PSR_ACC_CE=0`.

`provisional_ns_cleanup` is the previously reported
`rf_16_a8_e0p704546244` model. Its input pericentre is
2.3636300507767825 Rsun. At N-body time 0.001953125, CHAOS delivers
(a,e)=(7.833 Rsun, 0.69699) to EXPEL; the first COMENV call returns K=13,
e=0.001 and a=4.028 Rsun before the enforced second pass removes that
provisional state.

## Check the committed inputs

This check does not start NBODY6++GPU:

```bash
python3 tests/pulsar_registry_fixtures/verify_inputs.py
```

It verifies every checksum, input count, finite value, normalized BODY mass,
and the agreement between `dat.10`, `datsev.21`, and the namelist controls.

## Run a fixture

Build a serial executable compatible with the validation build:

```bash
./configure --with-par=1k --disable-gpu --disable-mpi --disable-hdf5 \
  --enable-simd=no --disable-omp
make -j
```

Run each case in its own empty directory. For example:

```bash
run_dir=$(mktemp -d)
cp tests/pulsar_registry_fixtures/ce_survive/{fixture.inp,dat.10,datsev.21} "$run_dir"/
cd "$run_dir"
/path/to/nbody6++ < fixture.inp > logfile 2> errfile
```

The detailed expected event sequence and final registry state for each case are
in `manifest.json`. This set covers the accepted COMENV/COAL, provisional-NS
cleanup, and CMBODY/MIX formation/migration paths. EXPEL2 chain-CE retry is
outside this fixture set.
