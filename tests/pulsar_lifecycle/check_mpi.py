#!/usr/bin/env python3
"""Check one-rank restart integration and two-rank rejection of bad dumps.

Requires matching serial/MPI dump ABI and results produced by run.py.
This is not a multi-rank integration regression (see README).
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess

from run import checkpoints, events, HERE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--serial-results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec')
    args = parser.parse_args()
    exe, base, out = args.executable.resolve(), args.serial_results.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    template = (HERE / 'rlof_ns_to_bh/restart.inp').read_text()
    failures = {
        'missing_record': 'missing pulsar dump record',
        'truncated_record': 'invalid pulsar dump state',
        'bad_count': 'invalid pulsar dump state',
        'stale_bh': 'inconsistent pulsar restart',
        'missing_ns': 'inconsistent pulsar restart',
        'duplicate': 'inconsistent pulsar restart',
        'disabled_history': 'cannot enable pulsars on a restart without PSR history',
        'invalid_restart_flags': 'pulsar evolution requires KZ(19)>=3',
    }
    reports = {}
    for label in ['rlof_ns_to_bh', 'post_collapse_restart'] + list(failures):
        directory = out / label; directory.mkdir()
        shutil.copy2(base / label / 'comm.1', directory / 'comm.1')
        inp = template
        ranks = 2 if label in failures else 1
        if label == 'post_collapse_restart':
            inp = inp.replace('TCRIT=0.2', 'TCRIT=0.02')
        if label == 'disabled_history':
            inp = inp.replace("Level='C'", "KZ(29)=2,Level='C'")
        if label == 'invalid_restart_flags':
            inp = inp.replace("Level='C'", "KZ(19)=2,Level='C'")
        with (directory / 'run.log').open('w') as log:
            result = subprocess.run(shlex.split(args.launcher) + ['-n', str(ranks), str(exe)],
                                    input=inp, text=True, cwd=directory, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=45,
                                    env=dict(os.environ, OMP_NUM_THREADS='1'))
        log = (directory / 'run.log').read_text(errors='replace')
        if label in failures:
            # MPI implementations differ in their launcher diagnostics.
            assert result.returncode != 0 and failures[label] in log, label
        else:
            assert result.returncode == 0 and 'END RUN' in log and 'FATAL' not in log, label
            final = checkpoints(directory)[1][-1]
            assert final['count'] == 0 and final['live'][1][0] == 14, label
            assert sum(row[1] == '6' for row in events(directory)) == int(label == 'rlof_ns_to_bh')
        reports[label] = dict(exit=result.returncode, ranks=ranks, status='PASS')
        print(f'{label}: MPI {ranks} rank(s) PASS', flush=True)
    (out / 'report.json').write_text(json.dumps(reports, indent=2) + '\n')


if __name__ == '__main__':
    main()
