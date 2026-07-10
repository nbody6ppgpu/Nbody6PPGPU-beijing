#!/usr/bin/env python3
"""Standard-library tests for the benchmark command-line tool."""

import importlib.util
import os
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name('benchmark.py')
SPEC = importlib.util.spec_from_file_location('nbody_benchmark', MODULE_PATH)
benchmark = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(benchmark)


class BenchmarkTests(unittest.TestCase):
    def test_parameter_combinations_and_directory_names(self):
        config = {
            'particle_number': '1k,2k',
            'node': '1',
            'mpi_per_node': '1,2',
            'gpu_per_node': '0',
            'openmp_thread_per_mpi': '4',
            'nbody_time': '1',
        }
        combinations = benchmark.generate_parameter_combinations(config)
        self.assertEqual(4, len(combinations))
        varying = benchmark.detect_varying_parameters(combinations)
        self.assertEqual(
            'N1k-1mpi', benchmark.generate_directory_name(combinations[0], varying)
        )
        self.assertEqual(
            'N2k-2mpi', benchmark.generate_directory_name(combinations[-1], varying)
        )

    def test_timing_result_parser(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / 'N1k-2mpi-4omp'
            run_dir.mkdir()
            output = run_dir / 'run.out'
            timing_values = ['1.0'] * 49
            output.write_text(
                'ADJUST TIME 2.0\n'
                + '0 1 1000 '
                + ' '.join(timing_values)
                + ' 0 0\n',
                encoding='utf-8',
            )
            rows = benchmark.extract_time_from_out_file(output)
            self.assertEqual(1, len(rows))
            self.assertEqual(1000, rows[0]['particle_number'])
            self.assertEqual(2, rows[0]['mpi_per_node'])
            self.assertEqual(4, rows[0]['openmp_thread_per_mpi'])
            self.assertEqual(2.0, rows[0]['NBTime'])

    def test_nonzero_exit_is_failure_even_with_completion_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            executable = run_dir / 'fake-nbody'
            executable.write_text('#!/bin/sh\necho "END RUN"\nexit 7\n', encoding='utf-8')
            executable.chmod(0o755)
            (run_dir / 'input.inp').write_text('', encoding='utf-8')

            success, output = benchmark.run_simulation_local(
                run_dir, executable, 'input.inp', 1, 2, False
            )
            self.assertFalse(success)
            self.assertTrue(Path(output).exists())


if __name__ == '__main__':
    unittest.main()
