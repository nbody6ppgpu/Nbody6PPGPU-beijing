#!/usr/bin/env python3
"""Serial lifecycle regressions; requires a completed 1k, no-SIMD build.

The generated restart is synthetic test state, not a physical birth fixture.
No production executable is instrumented; the normal integrator evolves it.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def records(path):
    result = []
    with path.open('rb') as f:
        while marker := f.read(4):
            size, = struct.unpack('<i', marker)
            assert size >= 0, path
            payload = f.read(size)
            assert len(payload) == size and f.read(4) == marker, path
            result.append(bytearray(payload))
    return result


def write_records(path, data):
    with path.open('wb') as f:
        for payload in data:
            marker = struct.pack('<i', len(payload))
            f.write(marker + payload + marker)


def dump(path):
    """This reader is specific to this tree's little-endian serial ABI."""
    rec = records(path)
    assert len(rec) == 10, (path, len(rec))
    ntot, pairs, _, n = struct.unpack_from('<4i', rec[1])
    params = struct.unpack_from('<265d', rec[1], 1024)
    particles = rec[2]
    assert len(particles) == 564 * ntot, path
    types = struct.unpack_from(f'<{ntot}i', particles, 552 * ntot)
    names = struct.unpack_from(f'<{ntot}i', particles, 560 * ntot)
    masses = struct.unpack_from(f'<{ntot}d', particles, 144 * ntot)
    live = {names[i]: (types[i], masses[i] * params[7])
            for i in range(n) if masses[i] > 0 and names[i] > 0}
    count, = struct.unpack('<i', rec[-2])
    assert len(rec[-1]) == 76 * count, path
    ints = struct.unpack_from(f'<{3 * count}i', rec[-1])
    vals = struct.unpack_from(f'<{8 * count}d', rec[-1], 12 * count)
    assert all(math.isfinite(v) for v in vals), path
    ids = list(ints[:count])
    assert len(set(ids)) == count, path
    ns = {name: mass for name, (kind, mass) in live.items() if kind == 13}
    assert set(ids) == set(ns), (path, ids, ns)
    for i, name in enumerate(ids):
        assert ints[count + i] == 13 and abs(vals[i] - ns[name]) < 1e-8, path
    return dict(time=params[157], count=count, names=ids,
                masses=list(vals[:count]), live=live, pairs=pairs)


def checkpoints(directory):
    paths = sorted(directory.glob('comm.2_*'),
                   key=lambda p: float(p.name.split('_', 1)[1]))
    assert paths, directory
    return paths, [dump(p) for p in paths]


def events(directory):
    return [line.split() for line in (directory / 'pulsar.204').read_text().splitlines()
            if line.startswith('PULSAR')]


def run(exe, directory, inp, failure=None):
    with (directory / 'run.log').open('w') as log:
        proc = subprocess.run([str(exe)], cwd=directory, input=inp, text=True,
                              stdout=log, stderr=subprocess.STDOUT, timeout=45,
                              env=dict(os.environ, OMP_NUM_THREADS='1'))
    log = (directory / 'run.log').read_text(errors='replace')
    if failure:
        assert proc.returncode != 0 and failure in log, (directory, proc.returncode)
    else:
        assert proc.returncode == 0, directory
        assert not any(x in log for x in ('FATAL ERROR', 'Fortran runtime error')), directory
    return log


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, default=ROOT / 'build/nbody6++')
    parser.add_argument('--output', type=Path, required=True, help='New directory only')
    parser.add_argument('--fc', default='gfortran', help='Compiler matching the serial build')
    args = parser.parse_args()
    exe, out = args.executable.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    common = [args.fc, '-O2', '-fPIC', '-mcmodel=large',
              '-I' + str(ROOT / 'extra_inc/nompi'), '-I' + str(ROOT / 'include')]
    objects = sorted((ROOT / 'build').glob('*.o'))
    assert objects and not list((ROOT / 'build').glob('*_mpi.o')), 'Use a clean serial build'
    for name, source in [('test_lifecycle', HERE / 'test_lifecycle.f'),
                         ('test_pulsar_invariants', ROOT / 'tests/pulsar_physics/test_pulsar_invariants.f')]:
        subprocess.run(common + [str(source), str(ROOT / 'build/pulsar.o'),
                                 str(ROOT / 'build/ran2.o'), '-o', str(out / name)], check=True)
        run(out / name, out, '')
        shutil.move(out / 'run.log', out / (name + '.log'))
        print(name + ': PASS', flush=True)
    seed_exe = out / 'seed_restart'
    subprocess.run(common + [str(HERE / 'seed_restart.f')] +
                   [str(p) for p in objects if p.name != 'nbody6.o'] +
                   ['-lstdc++', '-o', str(seed_exe)], check=True)
    fixture = HERE / 'rlof_ns_to_bh'
    seed = out / 'seed'; seed.mkdir()
    for name in ('dat.10', 'datsev.21'):
        shutil.copy2(fixture / name, seed / name)
    seed_input = (fixture / 'seed.inp').read_text()
    run(seed_exe, seed, seed_input)
    saved = seed / 'comm.1_0.000'
    assert dump(saved)['names'] == [1]
    restart_input = (fixture / 'restart.inp').read_text()

    def restart(label, data, inp=restart_input, failure=None):
        directory = out / label; directory.mkdir()
        write_records(directory / 'comm.1', data)
        log = run(exe, directory, inp, failure)
        if not failure:
            assert 'END RUN' in log, directory
        return directory, log

    collapse, log = restart('rlof_ns_to_bh', records(saved))
    paths, states = checkpoints(collapse)
    assert 'NEW ROCHE' in log and states[0]['count'] == 1
    final = states[-1]
    assert final['count'] == 0 and final['live'][1][0] == 14
    assert final['live'][1][1] > 2.5 and final['time'] >= .2 - 1e-12
    bank = events(collapse)
    retired = [r for r in bank if r[1] == '6']
    assert len(retired) == 1 and retired[0][3] == '1', bank
    assert bank[-1] == retired[0] and not any(r[1] == '1' for r in bank)
    print('RLOF NS -> BH, one retirement, all checkpoints: PASS', flush=True)
    resumed, _ = restart('post_collapse_restart', records(paths[-1]),
                         restart_input.replace('TCRIT=0.2', 'TCRIT=0.02'))
    assert all(s['count'] == 0 and s['live'][1][0] == 14
               for s in checkpoints(resumed)[1])
    assert not events(resumed)
    print('Post-collapse restart, no recreated NS/birth: PASS', flush=True)

    # Corrupt only the relevant record; never infer or repair birth history.
    for label in ('missing_record', 'truncated_record', 'bad_count', 'stale_bh', 'missing_ns', 'duplicate', 'disabled_history'):
        data = records(saved)
        inp = restart_input
        failure = 'inconsistent pulsar restart'
        if label == 'missing_record':
            data = data[:-2]; failure = 'missing pulsar dump record'
        elif label == 'truncated_record':
            data[-1] = data[-1][:4]; failure = 'invalid pulsar dump state'
        elif label == 'bad_count':
            data[-2] = struct.pack('<i', -1); failure = 'invalid pulsar dump state'
        elif label == 'stale_bh':
            ntot, = struct.unpack_from('<i', data[1])
            struct.pack_into('<i', data[2], 552 * ntot, 14)
        elif label == 'missing_ns':
            data[-2], data[-1] = struct.pack('<i', 0), b''
        elif label == 'duplicate':
            ints = struct.unpack_from('<3i', data[-1])
            vals = struct.unpack_from('<8d', data[-1], 12)
            data[-2] = struct.pack('<i', 2)
            data[-1] = struct.pack('<6i', *(v for v in ints for _ in range(2)))
            data[-1] += struct.pack('<16d', *(v for v in vals for _ in range(2)))
        else:
            struct.pack_into('<i', data[1], 228, 0)  # IA(55) = KZ(29)
            data = data[:-2]
            inp = inp.replace("Level='C'", "KZ(29)=2,Level='C'")
            failure = 'cannot enable pulsars on a restart without PSR history'
        restart(label, data, inp, failure)
    restart('invalid_restart_flags', records(saved),
            restart_input.replace("Level='C'", "KZ(19)=2,Level='C'"),
            'pulsar evolution requires KZ(19)>=3')
    print('Eight invalid restart cases rejected: PASS', flush=True)

    # Move INPULSAR into the normal fresh-start namelist position.
    fresh, control = seed_input.split('&INPULSAR', 1)
    fresh = fresh.replace('0 2\nKZ(31:', '2 2\nKZ(31:')
    fresh = fresh.replace('&INDATA', '&INPULSAR' + control + '\n&INDATA')
    for label, inp, failure in [
        ('initial_ns_rejected', fresh, 'initial NS pulsar state unsupported'),
        ('invalid_fresh_flags', fresh.replace('4 6\nKZ(21:', '2 6\nKZ(21:'),
         'pulsar evolution requires KZ(19)>=3'),
        ('disabled_initial_ns', seed_input, None),
    ]:
        directory = out / label; directory.mkdir()
        for name in ('dat.10', 'datsev.21'):
            shutil.copy2(fixture / name, directory / name)
        log = run(exe, directory, inp, failure)
        if not failure:
            assert 'END RUN' in log
    print('Initial NS and fresh flag policies, disabled-mode compatibility: PASS', flush=True)

    # Existing fixtures: the calibrated provisional case is only a smoke test.
    fixtures = ROOT / 'tests/pulsar_registry_fixtures'
    manifest = json.loads((fixtures / 'manifest.json').read_text())
    report = dict(executable=str(exe), sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                  compiler=subprocess.check_output([args.fc, '--version'], text=True).splitlines()[0],
                  collapse=final, fixtures={})
    for label, spec in manifest['fixtures'].items():
        directory = out / label; directory.mkdir()
        for name in ('dat.10', 'datsev.21'):
            shutil.copy2(fixtures / label / name, directory / name)
        log = run(exe, directory, (fixtures / label / 'fixture.inp').read_text())
        assert 'END RUN' in log
        paths, states = checkpoints(directory)
        final = states[-1]
        births = sum(r[1] == '1' for r in events(directory))
        assert not list(directory.glob('fort.2*'))
        if not spec.get('environment_sensitive'):
            expected = spec['expected']
            assert final['count'] == expected['final_ns_count']
            assert final['names'] == [expected['final_ns_name']]
            assert abs(final['masses'][0] - expected['final_ns_mass_msun']) < 1e-8
            assert births == 1
            if label.startswith('ce_'):
                assert 'ENFORCED CE' in log and 'BINARY CE :' in log
                assert ('BINARY COAL:' in log) == (label == 'ce_merge')
            else:
                rows = [r.split() for r in (directory / 'coll.13').read_text().splitlines()]
                assert any(r[3:6] == ['12', '11', '13'] for r in rows)
                if label == 'existing_ns_collision_migration':
                    assert any(r[3:6] == ['13', '1', '13'] for r in rows)
                    moved = [r[3] for r in events(directory) if r[1] == '5']
                    assert moved == ['1', '3']
            status = 'PASS'
        else:
            status = 'SMOKE ONLY (environment-sensitive target)'
        report['fixtures'][label] = dict(status=status, final=final, checkpoints=len(states), births=births)
        print(label + ': ' + status, flush=True)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Results: ' + str(out), flush=True)


if __name__ == '__main__':
    main()
