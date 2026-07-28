#!/usr/bin/env python3
"""Regression test for cumulative KZ(7)>=4 Lagrangian statistics."""

from __future__ import annotations

import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


NLAGR = 18
ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "examples" / "input_files" / "N1k_1Myr.inp"


def make_input(kz7: int) -> str:
    source = INPUT.read_text(encoding="utf-8")

    def replace_kz(match: re.Match[str]) -> str:
        values = match.group(1).split()
        if len(values) != 10:
            raise AssertionError("Unexpected KZ(1:10) input layout")
        values[6] = str(kz7)
        return "KZ(1:10)= " + " ".join(values)

    source, count = re.subn(
        r"KZ\(1:10\)=\s*([^\n]+)", replace_kz, source, count=1
    )
    if count != 1:
        raise AssertionError("KZ(1:10) was not found")
    source, count = re.subn(r"TCRIT=100\.0", "TCRIT=0.01", source, count=1)
    if count != 1:
        raise AssertionError("TCRIT was not found")
    return source


def run_case(binary: Path, kz7: int, parent: Path) -> tuple[list[float], str]:
    run_dir = parent / f"kz{kz7}"
    run_dir.mkdir()
    environment = os.environ.copy()
    environment.update(
        OMP_NUM_THREADS="4",
        OMP_STACKSIZE="512M",
    )
    result = subprocess.run(
        [str(binary)],
        input=make_input(kz7),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=run_dir,
        env=environment,
        timeout=120,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"KZ(7)={kz7} run failed ({result.returncode}):\n"
            + result.stdout[-4000:]
        )
    if "END RUN" not in result.stdout:
        raise AssertionError(f"KZ(7)={kz7} did not reach END RUN")

    records = [
        line
        for line in (run_dir / "lagr.7").read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#") and not line.startswith("TIME")
    ]
    if not records:
        raise AssertionError(f"KZ(7)={kz7} produced no lagr.7 data")
    values = [float(value.replace("D", "E")) for value in records[0].split()]
    if len(values) != 284:
        raise AssertionError(f"Expected 284 lagr.7 columns, got {len(values)}")
    return values, result.stdout


def group(record: list[float], start: int) -> list[float]:
    return record[start : start + NLAGR]


def validate_disabled_binary_radii(record: list[float]) -> None:
    for name, values in {
        "RSLAGR": group(record, 20),
        "RBLAGR": group(record, 38),
    }.items():
        if any(value != 0.0 for value in values):
            raise AssertionError(f"KZ(8)=0 {name} columns must be zero")


def assert_close(name: str, actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=2.0e-10, abs_tol=2.0e-12):
        raise AssertionError(
            f"{name}: actual={actual:.17e}, expected={expected:.17e}"
        )


def validate_main_output(output: str) -> None:
    for label in ("SIGR2:", "SIGT2:", "VROT:"):
        line = next((line for line in output.splitlines() if label in line), None)
        if line is None:
            raise AssertionError(f"Main output is missing {label}")
        values = [
            float(value.replace("D", "E"))
            for value in line.split(label, maxsplit=1)[1].split()
        ]
        if len(values) != NLAGR + 1:
            raise AssertionError(f"Main output {label} has {len(values)} values")
        if not all(math.isfinite(value) and abs(value) < 10.0 for value in values):
            raise AssertionError(f"Main output {label} is non-finite or implausible")


def validate_cumulative(shell: list[float], cumulative: list[float]) -> None:
    # Column offsets follow the header written by lagr.f.
    shell_mean_mass = group(shell, 56)
    shell_count = [int(value) for value in group(shell, 75)]
    cumulative_mean_mass = group(cumulative, 56)
    cumulative_count = [int(value) for value in group(cumulative, 75)]

    shell_fields = {
        "vx": group(shell, 94),
        "vy": group(shell, 113),
        "vz": group(shell, 132),
        "vr": group(shell, 170),
        "vrot": group(shell, 265),
    }
    cumulative_fields = {
        "vx": group(cumulative, 94),
        "vy": group(cumulative, 113),
        "vz": group(cumulative, 132),
        "vr": group(cumulative, 170),
        "vrot": group(cumulative, 265),
    }
    shell_sig2 = group(shell, 208)
    shell_sigr2 = group(shell, 227)
    shell_sigt2 = group(shell, 246)
    shell_vt = group(shell, 189)
    cumulative_sig2 = group(cumulative, 208)
    cumulative_sigr2 = group(cumulative, 227)
    cumulative_sigt2 = group(cumulative, 246)
    cumulative_vt = group(cumulative, 189)
    cumulative_v = group(cumulative, 151)

    checked_groups = (
        shell_mean_mass,
        shell_sig2,
        shell_sigr2,
        shell_sigt2,
        shell_vt,
        cumulative_mean_mass,
        cumulative_sig2,
        cumulative_sigr2,
        cumulative_sigt2,
        cumulative_vt,
        cumulative_v,
        *shell_fields.values(),
        *cumulative_fields.values(),
    )
    if not all(math.isfinite(value) for values in checked_groups for value in values):
        raise AssertionError("Lagrangian statistics contain non-finite values")

    running_mass = 0.0
    running_count = 0
    weighted = {name: 0.0 for name in shell_fields}
    raw_sig2 = 0.0
    raw_sigr2 = 0.0
    raw_sigt2 = 0.0

    for index in range(NLAGR):
        mass = shell_mean_mass[index] * shell_count[index]
        running_mass += mass
        running_count += shell_count[index]
        if running_mass <= 0.0:
            continue

        for name, values in shell_fields.items():
            weighted[name] += mass * values[index]
            assert_close(
                f"{name}[{index}]",
                cumulative_fields[name][index],
                weighted[name] / running_mass,
            )

        expected_count = running_count
        if cumulative_count[index] != expected_count:
            raise AssertionError(
                f"count[{index}]: {cumulative_count[index]} != {expected_count}"
            )
        assert_close(
            f"mean_mass[{index}]",
            cumulative_mean_mass[index],
            running_mass / running_count,
        )

        mean_speed2 = sum(
            cumulative_fields[name][index] ** 2 for name in ("vx", "vy", "vz")
        )
        assert_close(f"speed[{index}]", cumulative_v[index], math.sqrt(mean_speed2))

        shell_speed2 = sum(
            shell_fields[name][index] ** 2 for name in ("vx", "vy", "vz")
        )
        raw_sig2 += mass * (3.0 * shell_sig2[index] + shell_speed2)
        expected_sig2 = (raw_sig2 / running_mass - mean_speed2) / 3.0
        assert_close(f"sig2[{index}]", cumulative_sig2[index], expected_sig2)

        raw_sigr2 += mass * (
            shell_sigr2[index] + shell_fields["vr"][index] ** 2
        )
        expected_sigr2 = (
            raw_sigr2 / running_mass - cumulative_fields["vr"][index] ** 2
        )
        assert_close(f"sigr2[{index}]", cumulative_sigr2[index], expected_sigr2)

        raw_sigt2 += mass * (2.0 * shell_sigt2[index] + shell_vt[index] ** 2)
        expected_sigt2 = (
            raw_sigt2 / running_mass - cumulative_vt[index] ** 2
        ) / 2.0
        assert_close(f"sigt2[{index}]", cumulative_sigt2[index], expected_sigt2)

    for name, values in {
        "SIGR2": cumulative_sigr2,
        "SIGT2": cumulative_sigt2,
        "VROT": cumulative_fields["vrot"],
    }.items():
        if not all(math.isfinite(value) and abs(value) < 10.0 for value in values):
            raise AssertionError(f"{name} is non-finite or implausible")


def main() -> None:
    binary = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else ROOT / "build" / "nbody6++.sse"
    )
    if not binary.is_file():
        raise SystemExit(f"Binary not found: {binary}")
    with tempfile.TemporaryDirectory(prefix="nbody-lagr-") as temporary:
        shell, _ = run_case(binary, 3, Path(temporary))
        validate_disabled_binary_radii(shell)
        for kz7 in (4, 5):
            cumulative, output = run_case(binary, kz7, Path(temporary))
            validate_disabled_binary_radii(cumulative)
            validate_cumulative(shell, cumulative)
            validate_main_output(output)
    print("KZ(7)>=4 statistics and KZ(8)=0 radii passed")


if __name__ == "__main__":
    main()
