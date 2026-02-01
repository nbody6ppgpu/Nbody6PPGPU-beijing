#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Benchmark tool for NBODY6++GPU simulations.

This script enables:
1. Running benchmark simulations with various parameter combinations
2. Collecting and aggregating benchmark results

Author: NBODY6++GPU Team
"""

import argparse
import glob
import logging
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False


# Configure logging
logger = logging.getLogger('benchmark')
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def parse_particle_number(value: str) -> int:
    """
    Parse particle number with optional k/m suffix.
    
    Args:
        value: String like "50k", "100k", "0.5m", "2m", or plain integer
        
    Returns:
        Integer particle count
    """
    value = value.strip().lower()
    if value.endswith('m'):
        return int(float(value[:-1]) * 1_000_000)
    elif value.endswith('k'):
        return int(float(value[:-1]) * 1_000)
    else:
        return int(float(value))


def parse_comma_separated_ints(value: str) -> List[int]:
    """Parse comma-separated integers."""
    if not value:
        return []
    return [int(x.strip()) for x in value.split(',')]


def parse_comma_separated_particles(value: str) -> List[int]:
    """Parse comma-separated particle numbers with k/m suffixes."""
    if not value:
        return []
    return [parse_particle_number(x.strip()) for x in value.split(',')]


def parse_openmp_threads(value: str, mpi_per_node: int) -> List[Union[int, str]]:
    """
    Parse OpenMP thread specification.
    
    Args:
        value: String like "2,4" or "max"
        mpi_per_node: MPI processes per node (used for 'max' calculation)
        
    Returns:
        List of thread counts or ['max']
    """
    value = value.strip().lower()
    if value == 'max':
        return ['max']
    return [int(x.strip()) for x in value.split(',')]


def get_physical_cores() -> int:
    """Get the number of physical CPU cores."""
    try:
        # Try to get physical cores (not hyperthreaded)
        with open('/proc/cpuinfo', 'r') as f:
            content = f.read()
        # Count unique physical CPU cores
        physical_ids = set()
        core_ids = set()
        current_physical = None
        for line in content.split('\n'):
            if line.startswith('physical id'):
                current_physical = line.split(':')[1].strip()
            elif line.startswith('core id') and current_physical is not None:
                core_id = line.split(':')[1].strip()
                physical_ids.add((current_physical, core_id))
        if physical_ids:
            return len(physical_ids)
    except Exception:
        pass
    # Fallback to os.cpu_count()
    return os.cpu_count() or 1


def calculate_openmp_threads(mpi_per_node: int) -> int:
    """Calculate OpenMP threads when 'max' is specified."""
    physical_cores = get_physical_cores()
    return physical_cores // mpi_per_node


def format_particle_number(n: int) -> str:
    """Format particle number for directory naming."""
    if n >= 1_000_000 and n % 1_000_000 == 0:
        return f"{n // 1_000_000}m"
    elif n >= 1_000 and n % 1_000 == 0:
        return f"{n // 1_000}k"
    elif n >= 1_000_000:
        return f"{n / 1_000_000:.2f}m".rstrip('0').rstrip('.')
    elif n >= 1_000:
        return f"{n / 1_000:.2f}k".rstrip('0').rstrip('.')
    else:
        return str(n)


def generate_example_params() -> Dict[str, Any]:
    """Generate example benchmark parameters."""
    return {
        'particle_number': '50k,100k,0.5m,2m',
        'node': '1,2',
        'mpi_per_node': '1,2,4',
        'gpu_per_node': '4',
        'openmp_thread_per_mpi': '2,4',
        'nbody_time': '1',
        'disable_mpi': False,
        'expert': False,
    }


def load_yaml_config(filepath: str) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    if not YAML_AVAILABLE:
        logger.error("PyYAML is not installed. Please install it with: pip install pyyaml")
        sys.exit(1)
    
    with open(filepath, 'r') as f:
        return yaml.safe_load(f)


def save_yaml_config(filepath: str, config: Dict[str, Any]) -> None:
    """Save configuration to YAML file."""
    if not YAML_AVAILABLE:
        logger.error("PyYAML is not installed. Please install it with: pip install pyyaml")
        sys.exit(1)
    
    with open(filepath, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)


def find_executable(code_path: Path) -> Optional[Path]:
    """
    Find the NBODY6++ executable.
    
    Args:
        code_path: Path to the code repository
        
    Returns:
        Path to executable or None if not found
    """
    build_dir = code_path / 'build'
    pattern = str(build_dir / 'nbody6++.*')
    executables = sorted(glob.glob(pattern))
    if executables:
        return Path(executables[0])
    return None


def modify_input_file(
    input_path: Path,
    output_path: Path,
    n: int,
    tcrit: float
) -> None:
    """
    Modify input file with new parameters.
    
    Args:
        input_path: Path to base input file
        output_path: Path to write modified file
        n: Number of particles
        tcrit: Simulation time (TCRIT)
    """
    with open(input_path, 'r') as f:
        content = f.read()
    
    # Calculate DTMIN and RMIN
    # DTMIN * N = 0.6 => DTMIN = 0.6 / N
    # RMIN * N = 10.0 => RMIN = 10.0 / N
    dtmin = 0.6 / n
    rmin = 10.0 / n
    
    # Format in Fortran-compatible scientific notation
    def fortran_exp(val: float) -> str:
        """Format number in Fortran scientific notation (e.g., 6.0E-07)."""
        if val == 0:
            return "0.0E+00"
        exp = 0
        mantissa = val
        if abs(mantissa) >= 10:
            while abs(mantissa) >= 10:
                mantissa /= 10
                exp += 1
        elif abs(mantissa) < 1 and mantissa != 0:
            while abs(mantissa) < 1:
                mantissa *= 10
                exp -= 1
        sign = '+' if exp >= 0 else '-'
        return f"{mantissa:.1f}E{sign}{abs(exp):02d}"
    
    dtmin_str = fortran_exp(dtmin)
    rmin_str = fortran_exp(rmin)
    tcrit_str = f"{tcrit:.1f}"
    
    # Replace N= value (must be at start of line or after comma/space, not part of NRUN etc.)
    content = re.sub(r'(?<![A-Za-z])N=\d+', f'N={n}', content)
    
    # Replace TCRIT= value
    content = re.sub(r'(?<![A-Za-z])TCRIT=[\d.E+\-]+', f'TCRIT={tcrit_str}', content, flags=re.IGNORECASE)
    
    # Replace DTMIN= value
    content = re.sub(r'(?<![A-Za-z])DTMIN=[\d.E+\-]+', f'DTMIN={dtmin_str}', content, flags=re.IGNORECASE)
    
    # Replace RMIN= value (must not match GMIN, DTMIN etc.)
    content = re.sub(r'(?<![A-Za-z])RMIN=[\d.E+\-]+', f'RMIN={rmin_str}', content, flags=re.IGNORECASE)
    
    with open(output_path, 'w') as f:
        f.write(content)
    
    logger.debug(f"Modified input file: N={n}, TCRIT={tcrit_str}, DTMIN={dtmin_str}, RMIN={rmin_str}")


def modify_sbatch_file(
    sbatch_path: Path,
    output_path: Path,
    nodes: int,
    gpu_per_node: int,
    mpi_per_node: int,
    openmp_threads: int,
    exec_path: Path
) -> None:
    """
    Modify sbatch file with new parameters.
    
    Args:
        sbatch_path: Path to base sbatch file
        output_path: Path to write modified file
        nodes: Number of nodes
        gpu_per_node: GPUs per node
        mpi_per_node: MPI processes per node
        openmp_threads: OpenMP threads per MPI process
        exec_path: Path to executable
    """
    with open(sbatch_path, 'r') as f:
        content = f.read()
    
    # Replace SBATCH parameters
    content = re.sub(r'#SBATCH --nodes=\d+', f'#SBATCH --nodes={nodes}', content)
    content = re.sub(r'#SBATCH --gres=gpu:\d+', f'#SBATCH --gres=gpu:{gpu_per_node}', content)
    content = re.sub(r'#SBATCH --ntasks-per-node=\d+', f'#SBATCH --ntasks-per-node={mpi_per_node}', content)
    content = re.sub(r'#SBATCH --cpus-per-task=\d+', f'#SBATCH --cpus-per-task={openmp_threads}', content)
    
    # Replace EXEC= line
    content = re.sub(r'EXEC=.*', f'EXEC={exec_path}', content)
    
    with open(output_path, 'w') as f:
        f.write(content)
    
    logger.debug(f"Modified sbatch file: nodes={nodes}, gpus={gpu_per_node}, mpi={mpi_per_node}, omp={openmp_threads}")


def generate_parameter_combinations(config: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Generate all parameter combinations for benchmarking.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        List of parameter combination dictionaries
    """
    particle_numbers = parse_comma_separated_particles(config.get('particle_number', ''))
    nodes = parse_comma_separated_ints(str(config.get('node', '1')))
    mpi_per_node_list = parse_comma_separated_ints(str(config.get('mpi_per_node', '1')))
    gpu_per_node_list = parse_comma_separated_ints(str(config.get('gpu_per_node', '1')))
    nbody_times = parse_comma_separated_ints(str(config.get('nbody_time', '1')))
    
    omp_spec = str(config.get('openmp_thread_per_mpi', '1'))
    
    combinations = []
    
    for n in particle_numbers:
        for node in nodes:
            for mpi in mpi_per_node_list:
                for gpu in gpu_per_node_list:
                    for nbtime in nbody_times:
                        # Handle OpenMP threads
                        if omp_spec.strip().lower() == 'max':
                            omp_threads = calculate_openmp_threads(mpi)
                            omp_values = [omp_threads]
                        else:
                            omp_values = parse_comma_separated_ints(omp_spec)
                        
                        for omp in omp_values:
                            combinations.append({
                                'particle_number': n,
                                'node': node,
                                'mpi_per_node': mpi,
                                'gpu_per_node': gpu,
                                'openmp_thread_per_mpi': omp,
                                'nbody_time': nbtime,
                            })
    
    return combinations


def generate_directory_name(
    params: Dict[str, Any],
    varying_params: Dict[str, bool]
) -> str:
    """
    Generate directory name based on varying parameters.
    
    Args:
        params: Parameter values for this combination
        varying_params: Which parameters have multiple values
        
    Returns:
        Directory name string
    """
    n = params['particle_number']
    parts = [f"N{format_particle_number(n)}"]
    
    # Always include node if varying
    if varying_params.get('node', False):
        parts.append(f"{params['node']}node")
    
    if varying_params.get('mpi_per_node', False):
        parts.append(f"{params['mpi_per_node']}mpi")
    
    if varying_params.get('gpu_per_node', False):
        parts.append(f"{params['gpu_per_node']}gpu")
    
    if varying_params.get('openmp_thread_per_mpi', False):
        parts.append(f"{params['openmp_thread_per_mpi']}omp")
    
    if varying_params.get('nbody_time', False):
        parts.append(f"{params['nbody_time']}T")
    
    return '-'.join(parts)


def detect_varying_parameters(combinations: List[Dict[str, Any]]) -> Dict[str, bool]:
    """
    Detect which parameters vary across combinations.
    
    Args:
        combinations: List of parameter combinations
        
    Returns:
        Dictionary indicating which parameters vary
    """
    if not combinations:
        return {}
    
    varying = {}
    first = combinations[0]
    
    for key in first:
        values = set(combo[key] for combo in combinations)
        varying[key] = len(values) > 1
    
    return varying


def log_code_info(code_path: Path, run_dir: Path) -> None:
    """
    Log code information to run directory.
    
    Args:
        code_path: Path to code repository
        run_dir: Path to run directory
    """
    log_file = run_dir / 'benchmark_info.log'
    
    with open(log_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("NBODY6++GPU Benchmark Information\n")
        f.write(f"Timestamp: {datetime.now().isoformat()}\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Code commit id:\n")
        try:
            result = subprocess.run(
                ['git', 'log', '-1', '--pretty=format:%h %s'],
                cwd=code_path,
                capture_output=True,
                text=True
            )
            f.write(result.stdout + "\n\n")
        except Exception as e:
            f.write(f"Error getting git info: {e}\n\n")
        
        f.write("Code build date:\n")
        try:
            exec_path = find_executable(code_path)
            if exec_path:
                result = subprocess.run(
                    ['ls', '-lh', str(exec_path)],
                    capture_output=True,
                    text=True
                )
                f.write(result.stdout + "\n")
            else:
                f.write("No executable found\n")
        except Exception as e:
            f.write(f"Error getting build info: {e}\n")
    
    # Copy config.log if exists
    config_log = code_path / 'config.log'
    if config_log.exists():
        shutil.copy(config_log, run_dir / 'config.log')
    
    logger.info(f"Code information logged to {log_file}")


def run_simulation_slurm(run_dir: Path, sbatch_file: str) -> Tuple[bool, str]:
    """
    Submit simulation job via SLURM.
    
    Args:
        run_dir: Directory containing sbatch file
        sbatch_file: Name of sbatch file
        
    Returns:
        Tuple of (success, message)
    """
    try:
        result = subprocess.run(
            ['sbatch', sbatch_file],
            cwd=run_dir,
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            logger.info(f"Submitted job in {run_dir}: {result.stdout.strip()}")
            return True, result.stdout.strip()
        else:
            logger.error(f"Failed to submit job: {result.stderr}")
            return False, result.stderr
    except Exception as e:
        logger.error(f"Error submitting job: {e}")
        return False, str(e)


def run_simulation_local(
    run_dir: Path,
    exec_path: Path,
    input_file: str,
    mpi_procs: int,
    openmp_threads: int,
    use_mpi: bool
) -> Tuple[bool, str]:
    """
    Run simulation locally without SLURM.
    
    Args:
        run_dir: Directory to run simulation in
        exec_path: Path to executable
        input_file: Input file name
        mpi_procs: Total MPI processes
        openmp_threads: OpenMP threads per process
        use_mpi: Whether to use MPI
        
    Returns:
        Tuple of (success, output_file_path)
    """
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    basename = Path(input_file).stem
    out_file = f"{basename}.{timestamp}.out"
    err_file = f"{basename}.{timestamp}.err"
    
    # Set environment
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(openmp_threads)
    env['OMP_STACKSIZE'] = '10G'
    env['OMP_PROC_BIND'] = 'true'
    
    # Build command
    if use_mpi and '.mpi' in str(exec_path):
        cmd = ['mpirun', '-n', str(mpi_procs), '--bind-to', 'none', str(exec_path)]
    else:
        cmd = [str(exec_path)]
    
    logger.info(f"Running: {' '.join(cmd)} in {run_dir}")
    logger.info(f"OMP_NUM_THREADS={openmp_threads}")
    
    try:
        with open(run_dir / input_file, 'r') as stdin_file:
            with open(run_dir / out_file, 'w') as stdout_file:
                with open(run_dir / err_file, 'w') as stderr_file:
                    result = subprocess.run(
                        cmd,
                        cwd=run_dir,
                        stdin=stdin_file,
                        stdout=stdout_file,
                        stderr=stderr_file,
                        env=env
                    )
        
        # Check for successful completion
        with open(run_dir / out_file, 'r') as f:
            content = f.read()
            if re.search(r'END RUN|TERMINATE', content):
                logger.info(f"Simulation completed successfully: {out_file}")
                return True, str(run_dir / out_file)
            else:
                logger.warning(f"Simulation may have failed (no END RUN found): {out_file}")
                return False, str(run_dir / out_file)
                
    except Exception as e:
        logger.error(f"Error running simulation: {e}")
        return False, str(e)


def collect_benchmark_results(run_dir: Path) -> Optional['pd.DataFrame']:
    """
    Collect benchmark results from output files.
    
    This function parses the benchmark timing data from simulation output files,
    similar to the logic in get_profile_time_csv.
    
    Args:
        run_dir: Directory containing benchmark runs
        
    Returns:
        DataFrame with collected results, or None if pandas not available
    """
    if not PANDAS_AVAILABLE:
        logger.error("Pandas is not installed. Please install it with: pip install pandas")
        return None
    
    # Headers based on get_profile_time_csv
    columns = [
        'run_dir', 'NBTime', 'rank', 'PE', 'N', 'Total', 'Inti.', 'Intgrt', 'Reg.',
        'Irr.', 'Predall', 'Pred.', 'Init.B.', 'Mdot', 'Move', 'Comm.I.', 'Comm.R.',
        'Send.I.', 'Send.R.', 'KS', 'Adjust', 'OUT', 'Barr.', 'Barr.P.', 'Barr.I.',
        'Barr.R.', 'Reg.GPU.S', 'Reg.GPU.P', 'Comm.Adj.', 'Mdot.Fic.', 'Mdot.Fc.',
        'Mdot.Pot.', 'Mdot.EC.', 'Sort.B.', 'HighV', 'KS.Init.B', 'KS.Int.S',
        'KS.Int.P', 'KS.Comm.', 'KS.Barr.', 'KS.Move', 'KS.Cmb.', 'KS.Insert',
        'KS.Init.', 'KS.Term.', 'Hiar.', 'KS.UP', 'KS.TP', 'TIDES3', 'GRRAD',
        'xtsub1', 'xtsub2', 'xnirrf', 'xnpred', 'itides3', 'igrrad',
        'node', 'mpi_per_node', 'gpu_per_node', 'openmp_thread_per_mpi',
        'particle_number', 'nbody_time'
    ]
    
    results = []
    
    # Find all subdirectories
    for subdir in run_dir.iterdir():
        if not subdir.is_dir():
            continue
        if subdir.name.startswith('.'):
            continue
        
        # Find .out files in this directory
        out_files = list(subdir.glob('*.out'))
        
        for out_file in out_files:
            try:
                with open(out_file, 'r') as f:
                    content = f.read()
                
                # Look for ADJUST block with timing data
                # Pattern based on get_profile_time_csv
                adjust_matches = re.findall(
                    r'ADJUST.*?TIME\s+([\d.DE+-]+).*?\n.*?(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+).*?0\.00E\+00\s+0\.00E\+00',
                    content,
                    re.DOTALL
                )
                
                if not adjust_matches:
                    # Try simpler pattern for benchmark data
                    # Look for lines with timing information
                    time_match = re.search(r'TIME\s+([\d.DE+-]+)', content)
                    if time_match:
                        nb_time = time_match.group(1).replace('D', 'E')
                    else:
                        continue
                
                # Parse directory name for configuration
                dir_name = subdir.name
                config = parse_directory_name(dir_name)
                
                # Find the main timing line
                # Looking for pattern like: rank PE N Total ...
                timing_pattern = r'^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+(\d+)\s+(\d+)'
                
                for line in content.split('\n'):
                    match = re.match(timing_pattern, line.strip())
                    if match:
                        # Extract TIME value from nearby ADJUST line
                        nb_time_val = 0.0
                        time_search = re.search(r'TIME\s+([\d.DE+-]+)', content)
                        if time_search:
                            nb_time_val = float(time_search.group(1).replace('D', 'E'))
                        
                        row = {
                            'run_dir': dir_name,
                            'NBTime': nb_time_val,
                            'rank': int(match.group(1)),
                            'PE': int(match.group(2)),
                            'N': int(match.group(3)),
                            'Total': float(match.group(4)),
                            'Inti.': float(match.group(5)),
                            'Intgrt': float(match.group(6)),
                            'Reg.': float(match.group(7)),
                            'Irr.': float(match.group(8)),
                            'Predall': float(match.group(9)),
                            'Pred.': float(match.group(10)),
                            'Init.B.': float(match.group(11)),
                            'Mdot': float(match.group(12)),
                            'Move': float(match.group(13)),
                            'Comm.I.': float(match.group(14)),
                            'Comm.R.': float(match.group(15)),
                            'Send.I.': float(match.group(16)),
                            'Send.R.': float(match.group(17)),
                            'KS': float(match.group(18)),
                            'Adjust': float(match.group(19)),
                            'OUT': float(match.group(20)),
                            'Barr.': float(match.group(21)),
                            'Barr.P.': float(match.group(22)),
                            'Barr.I.': float(match.group(23)),
                            'Barr.R.': float(match.group(24)),
                            'Reg.GPU.S': float(match.group(25)),
                            'Reg.GPU.P': float(match.group(26)),
                            'Comm.Adj.': float(match.group(27)),
                            'Mdot.Fic.': float(match.group(28)),
                            'Mdot.Fc.': float(match.group(29)),
                            'Mdot.Pot.': float(match.group(30)),
                            'Mdot.EC.': float(match.group(31)),
                            'Sort.B.': float(match.group(32)),
                            'HighV': float(match.group(33)),
                            'KS.Init.B': float(match.group(34)),
                            'KS.Int.S': float(match.group(35)),
                            'KS.Int.P': float(match.group(36)),
                            'KS.Comm.': float(match.group(37)),
                            'KS.Barr.': float(match.group(38)),
                            'KS.Move': float(match.group(39)),
                            'KS.Cmb.': float(match.group(40)),
                            'KS.Insert': float(match.group(41)),
                            'KS.Init.': float(match.group(42)),
                            'KS.Term.': float(match.group(43)),
                            'Hiar.': float(match.group(44)),
                            'KS.UP': float(match.group(45)),
                            'KS.TP': float(match.group(46)),
                            'TIDES3': float(match.group(47)),
                            'GRRAD': float(match.group(48)),
                            'xtsub1': float(match.group(49).replace('E', 'e')),
                            'xtsub2': float(match.group(50).replace('E', 'e')),
                            'xnirrf': float(match.group(51).replace('E', 'e')),
                            'xnpred': float(match.group(52).replace('E', 'e')),
                            'itides3': int(match.group(53)),
                            'igrrad': int(match.group(54)),
                        }
                        
                        # Add configuration from directory name
                        row.update(config)
                        results.append(row)
                        break
                        
            except Exception as e:
                logger.warning(f"Error processing {out_file}: {e}")
                continue
    
    if not results:
        logger.warning("No benchmark results found")
        return None
    
    df = pd.DataFrame(results)
    return df


def parse_directory_name(dir_name: str) -> Dict[str, Any]:
    """
    Parse configuration from directory name.
    
    Args:
        dir_name: Directory name like "N100k-2node-4mpi-2gpu-8omp"
        
    Returns:
        Configuration dictionary
    """
    config = {
        'node': 1,
        'mpi_per_node': 1,
        'gpu_per_node': 1,
        'openmp_thread_per_mpi': 1,
        'particle_number': 0,
        'nbody_time': 1,
    }
    
    # Parse N (particle number)
    n_match = re.search(r'N([\d.]+[km]?)', dir_name, re.IGNORECASE)
    if n_match:
        config['particle_number'] = parse_particle_number(n_match.group(1))
    
    # Parse node count
    node_match = re.search(r'(\d+)node', dir_name)
    if node_match:
        config['node'] = int(node_match.group(1))
    
    # Parse MPI per node
    mpi_match = re.search(r'(\d+)mpi', dir_name)
    if mpi_match:
        config['mpi_per_node'] = int(mpi_match.group(1))
    
    # Parse GPU per node
    gpu_match = re.search(r'(\d+)gpu', dir_name)
    if gpu_match:
        config['gpu_per_node'] = int(gpu_match.group(1))
    
    # Parse OpenMP threads
    omp_match = re.search(r'(\d+)omp', dir_name)
    if omp_match:
        config['openmp_thread_per_mpi'] = int(omp_match.group(1))
    
    # Parse nbody time
    time_match = re.search(r'(\d+)T', dir_name)
    if time_match:
        config['nbody_time'] = int(time_match.group(1))
    
    return config


def user_confirm(prompt: str) -> bool:
    """
    Ask user for confirmation.
    
    Args:
        prompt: Prompt message
        
    Returns:
        True if user confirms, False otherwise
    """
    while True:
        response = input(f"{prompt} [y/n]: ").strip().lower()
        if response in ('y', 'yes'):
            return True
        elif response in ('n', 'no'):
            return False
        print("Please enter 'y' or 'n'")


def validate_config(config: Dict[str, Any], slurm_mode: bool) -> bool:
    """
    Validate configuration.
    
    Args:
        config: Configuration dictionary
        slurm_mode: Whether running in SLURM mode
        
    Returns:
        True if valid, False otherwise
    """
    errors = []
    
    # Check required parameters
    if not config.get('particle_number'):
        errors.append("particle_number is required")
    
    if not config.get('exec_path'):
        errors.append("exec_path is required or must be auto-detected")
    elif not Path(config['exec_path']).exists():
        errors.append(f"Executable not found: {config['exec_path']}")
    
    if not config.get('input_base_path'):
        errors.append("input_base_path is required")
    elif not Path(config['input_base_path']).exists():
        errors.append(f"Input file not found: {config['input_base_path']}")
    
    if slurm_mode and not config.get('sbatch_base_path'):
        errors.append("sbatch_base_path is required for SLURM mode")
    elif slurm_mode and not Path(config['sbatch_base_path']).exists():
        errors.append(f"Sbatch file not found: {config['sbatch_base_path']}")
    
    # Check GPU warning
    if config.get('gpu_per_node') and not slurm_mode:
        logger.warning("gpu_per_node is specified but sbatch_base_path is not set. "
                      "GPU setting only takes effect in SLURM mode.")
    
    if errors:
        for error in errors:
            logger.error(error)
        return False
    
    return True


def setup_argument_parser() -> argparse.ArgumentParser:
    """Create and configure argument parser."""
    parser = argparse.ArgumentParser(
        description='Benchmark tool for NBODY6++GPU simulations.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate example parameter file
  %(prog)s --generate-example-param
  
  # Run benchmark with parameter file
  %(prog)s --param-file-path=benchmark_params.yaml
  
  # Run benchmark with command line arguments
  %(prog)s -N 50k,100k --node=1,2 --mpi-per-node=4
  
  # Collect results only
  %(prog)s --collect-only --param-file-path=benchmark_params.yaml
        """
    )
    
    # Parameter file options
    parser.add_argument(
        '--generate-example-param',
        action='store_true',
        help='Generate an example parameter file (YAML format)'
    )
    parser.add_argument(
        '--param-file-path',
        type=str,
        help='Path to benchmark parameter file (YAML format)'
    )
    
    # SLURM options
    parser.add_argument(
        '--sbatch-base-path',
        type=str,
        help='Path to base sbatch file. If not specified, runs without SLURM. '
             'Example file available at: code-path/examples/example.sbatch'
    )
    
    # Simulation parameters
    parser.add_argument(
        '-N', '--particle-number',
        type=str,
        help='Comma-separated particle numbers. Supports k/m suffix (e.g., 50k,100k,0.5m,2m)'
    )
    parser.add_argument(
        '--node',
        type=str,
        default='1',
        help='Comma-separated node counts (e.g., 1,2). Default: 1'
    )
    parser.add_argument(
        '--mpi-per-node',
        type=str,
        default='1',
        help='Comma-separated MPI processes per node (e.g., 1,2,4). Default: 1'
    )
    parser.add_argument(
        '--gpu-per-node',
        type=str,
        default='4',
        help='Comma-separated GPUs per node (e.g., 4). Only effective in SLURM mode. Default: 4'
    )
    parser.add_argument(
        '--openmp-thread-per-mpi',
        type=str,
        default='1',
        help='Comma-separated OpenMP threads per MPI, or "max" for auto-detection. Default: 1'
    )
    parser.add_argument(
        '--nbody-time',
        type=str,
        default='1',
        help='Comma-separated N-body simulation times (e.g., 1). Default: 1'
    )
    parser.add_argument(
        '--disable-mpi',
        action='store_true',
        help='Disable MPI execution (run single process)'
    )
    
    # Path options
    parser.add_argument(
        '--code-path',
        type=str,
        help='Path to code repository. Default: parent directory of this script'
    )
    parser.add_argument(
        '--input-base-path',
        type=str,
        help='Path to base input file. Default: code-path/examples/N1m_benchmark.inp'
    )
    parser.add_argument(
        '--exec-path',
        type=str,
        help='Path to executable. Auto-detected from code-path/build/ if not specified'
    )
    parser.add_argument(
        '--run-dir',
        type=str,
        help='Directory for benchmark runs. Default: code-path/benchmark_run'
    )
    
    # Execution options
    parser.add_argument(
        '--collect-only',
        action='store_true',
        help='Only collect results without running simulations'
    )
    parser.add_argument(
        '--expert',
        action='store_true',
        help='Skip confirmation prompts'
    )
    
    # Logging
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress non-error messages'
    )
    
    return parser


def merge_config(args: argparse.Namespace, yaml_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Merge command line arguments with YAML configuration.
    Command line arguments take precedence.
    
    Args:
        args: Parsed command line arguments
        yaml_config: Configuration from YAML file
        
    Returns:
        Merged configuration dictionary
    """
    config = yaml_config.copy()
    
    # Map argument names to config keys
    arg_mapping = {
        'sbatch_base_path': 'sbatch_base_path',
        'particle_number': 'particle_number',
        'node': 'node',
        'mpi_per_node': 'mpi_per_node',
        'gpu_per_node': 'gpu_per_node',
        'openmp_thread_per_mpi': 'openmp_thread_per_mpi',
        'nbody_time': 'nbody_time',
        'disable_mpi': 'disable_mpi',
        'code_path': 'code_path',
        'input_base_path': 'input_base_path',
        'exec_path': 'exec_path',
        'run_dir': 'run_dir',
        'expert': 'expert',
    }
    
    for arg_name, config_key in arg_mapping.items():
        arg_value = getattr(args, arg_name, None)
        if arg_value is not None:
            config[config_key] = arg_value
    
    return config


def main():
    """Main entry point."""
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # Configure logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    elif args.quiet:
        logger.setLevel(logging.ERROR)
    
    # Determine script directory and default paths
    script_dir = Path(__file__).resolve().parent
    default_code_path = script_dir.parent
    
    # Handle --generate-example-param
    if args.generate_example_param:
        example_params = generate_example_params()
        example_file = Path.cwd() / 'benchmark_params.yaml'
        save_yaml_config(str(example_file), example_params)
        logger.info(f"Generated example parameter file: {example_file}")
        return 0
    
    # Load YAML config if specified
    yaml_config = {}
    if args.param_file_path:
        logger.info(f"Loading parameters from: {args.param_file_path}")
        yaml_config = load_yaml_config(args.param_file_path)
    
    # Merge configurations
    config = merge_config(args, yaml_config)
    
    # Set default paths
    code_path = Path(config.get('code_path', default_code_path)).resolve()
    config['code_path'] = str(code_path)
    
    if not config.get('input_base_path'):
        config['input_base_path'] = str(code_path / 'examples' / 'N1m_benchmark.inp')
    
    if not config.get('run_dir'):
        config['run_dir'] = str(code_path / 'benchmark_run')
    
    if not config.get('exec_path'):
        exec_path = find_executable(code_path)
        if exec_path:
            config['exec_path'] = str(exec_path)
            logger.info(f"Auto-detected executable: {exec_path}")
        else:
            logger.error("Could not find executable. Please compile the code first with:")
            logger.error("  ./configure && make")
            logger.error("Or specify --exec-path manually")
            return 1
    
    # Determine if SLURM mode
    slurm_mode = bool(config.get('sbatch_base_path'))
    
    # Handle --collect-only
    if args.collect_only:
        run_dir = Path(config['run_dir'])
        if not run_dir.exists():
            logger.error(f"Run directory does not exist: {run_dir}")
            return 1
        
        logger.info("Collecting benchmark results...")
        df = collect_benchmark_results(run_dir)
        if df is not None:
            output_csv = run_dir / 'benchmark_results.csv'
            df.to_csv(output_csv, index=False)
            logger.info(f"Results saved to: {output_csv}")
            logger.info(f"Collected {len(df)} result entries")
        return 0
    
    # Validate configuration
    if not validate_config(config, slurm_mode):
        return 1
    
    # Create run directory
    run_dir = Path(config['run_dir'])
    run_dir.mkdir(parents=True, exist_ok=True)
    
    # Log code information
    log_code_info(code_path, run_dir)
    
    # Save current configuration
    config_file = run_dir / 'benchmark_config.yaml'
    save_yaml_config(str(config_file), config)
    logger.info(f"Configuration saved to: {config_file}")
    
    # Generate parameter combinations
    combinations = generate_parameter_combinations(config)
    if not combinations:
        logger.error("No parameter combinations generated. Check particle_number setting.")
        return 1
    
    logger.info(f"Generated {len(combinations)} parameter combinations")
    
    # Detect varying parameters for directory naming
    varying_params = detect_varying_parameters(combinations)
    
    # Prepare directories and files
    input_base_path = Path(config['input_base_path'])
    exec_path = Path(config['exec_path'])
    sbatch_base_path = Path(config['sbatch_base_path']) if slurm_mode else None
    
    prepared_dirs = []
    for params in combinations:
        dir_name = generate_directory_name(params, varying_params)
        sub_dir = run_dir / dir_name
        sub_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy and modify input file
        input_file = sub_dir / input_base_path.name
        modify_input_file(
            input_base_path,
            input_file,
            params['particle_number'],
            float(params['nbody_time'])
        )
        
        # Copy and modify sbatch file if SLURM mode
        if slurm_mode:
            sbatch_file = sub_dir / sbatch_base_path.name
            modify_sbatch_file(
                sbatch_base_path,
                sbatch_file,
                params['node'],
                params['gpu_per_node'],
                params['mpi_per_node'],
                params['openmp_thread_per_mpi'],
                exec_path
            )
        
        prepared_dirs.append({
            'dir': sub_dir,
            'params': params,
            'input_file': input_file.name,
            'sbatch_file': sbatch_base_path.name if slurm_mode else None,
        })
        
        logger.debug(f"Prepared: {dir_name}")
    
    logger.info(f"Prepared {len(prepared_dirs)} benchmark directories in {run_dir}")
    
    # Confirm before running
    if not config.get('expert', False):
        print("\nBenchmark configuration:")
        print(f"  Run directory: {run_dir}")
        print(f"  Executable: {exec_path}")
        print(f"  Parameter combinations: {len(combinations)}")
        print(f"  SLURM mode: {slurm_mode}")
        if not user_confirm("Start simulations?"):
            logger.info("Aborted by user")
            return 0
    
    # Run simulations
    results = []
    use_mpi = not config.get('disable_mpi', False)
    
    for prep in prepared_dirs:
        sub_dir = prep['dir']
        params = prep['params']
        
        if slurm_mode:
            # Submit SLURM job
            success, msg = run_simulation_slurm(sub_dir, prep['sbatch_file'])
            results.append({'dir': sub_dir, 'success': success, 'message': msg})
        else:
            # Run locally
            mpi_procs = params['node'] * params['mpi_per_node']
            success, out_file = run_simulation_local(
                sub_dir,
                exec_path,
                prep['input_file'],
                mpi_procs,
                params['openmp_thread_per_mpi'],
                use_mpi
            )
            results.append({'dir': sub_dir, 'success': success, 'output': out_file})
    
    # Report results
    success_count = sum(1 for r in results if r['success'])
    logger.info(f"\nSimulation results: {success_count}/{len(results)} successful")
    
    if not slurm_mode:
        # For local runs, check for failures
        failures = [r for r in results if not r['success']]
        if failures:
            logger.warning("Some simulations may have failed:")
            for f in failures:
                logger.warning(f"  {f['dir']}")
    
    # Collect results (for local runs)
    if not slurm_mode:
        if not config.get('expert', False):
            if not user_confirm("Collect benchmark results?"):
                logger.info("Skipping result collection")
                return 0
        
        logger.info("Collecting benchmark results...")
        df = collect_benchmark_results(run_dir)
        if df is not None:
            output_csv = run_dir / 'benchmark_results.csv'
            df.to_csv(output_csv, index=False)
            logger.info(f"Results saved to: {output_csv}")
            logger.info(f"Collected {len(df)} result entries")
    else:
        logger.info("SLURM jobs submitted. Run with --collect-only after jobs complete to gather results.")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
