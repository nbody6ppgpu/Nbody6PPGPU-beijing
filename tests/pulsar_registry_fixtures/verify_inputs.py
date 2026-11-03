#!/usr/bin/env python3
"""Verify the committed small-N fixture inputs without running a simulation."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REQUIRED_FILES = {"fixture.inp", "dat.10", "datsev.21"}


def fail(message: str) -> None:
    raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def namelist_number(text: str, name: str) -> float:
    pattern = rf"(?i)(?<![A-Z0-9_]){re.escape(name)}\s*=\s*([-+0-9.Ee]+)"
    match = re.search(pattern, text)
    if not match:
        fail(f"missing {name} in fixture.inp")
    return float(match.group(1))


def parse_numeric_rows(path: Path, columns: int) -> list[list[float]]:
    rows: list[list[float]] = []
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        if not raw.strip():
            continue
        fields = raw.split()
        if len(fields) != columns:
            fail(f"{path}: line {line_number} has {len(fields)} columns, expected {columns}")
        try:
            row = [float(field) for field in fields]
        except ValueError as exc:
            fail(f"{path}: line {line_number} is not numeric: {exc}")
        if not all(math.isfinite(value) for value in row):
            fail(f"{path}: line {line_number} contains a non-finite value")
        rows.append(row)
    return rows


def verify_case(case_name: str, case: dict[str, object]) -> None:
    case_dir = ROOT / case_name
    if not case_dir.is_dir():
        fail(f"missing fixture directory: {case_name}")

    actual_files = {
        path.name for path in case_dir.iterdir() if path.is_file()
    }
    if actual_files != REQUIRED_FILES:
        fail(f"{case_name}: files are {sorted(actual_files)}, expected {sorted(REQUIRED_FILES)}")

    expected_hashes = case["files"]
    if not isinstance(expected_hashes, dict):
        fail(f"{case_name}: invalid files entry in manifest")
    for filename in sorted(REQUIRED_FILES):
        actual = sha256(case_dir / filename)
        expected = expected_hashes.get(filename)
        if actual != expected:
            fail(f"{case_name}/{filename}: SHA256 {actual}, expected {expected}")

    inp = (case_dir / "fixture.inp").read_text()
    for name, expected in {
        "N": 50,
        "NFIX": 1,
        "NRAND": 43532,
        "KSTART": 1,
        "PSR_ACC_CE": 0,
    }.items():
        actual = namelist_number(inp, name)
        if actual != expected:
            fail(f"{case_name}: {name}={actual}, expected {expected}")

    for index, expected in {19: 4, 22: 2, 29: 2, 50: 1}.items():
        block_start = ((index - 1) // 10) * 10 + 1
        block_end = block_start + 9
        match = re.search(
            rf"KZ\({block_start}:{block_end}\)\s*=\s*([^\n]+)", inp, re.IGNORECASE
        )
        if not match:
            fail(f"{case_name}: missing KZ({block_start}:{block_end}) block")
        values = [int(value) for value in re.findall(r"[-+]?\d+", match.group(1))]
        if len(values) != 10:
            fail(f"{case_name}: malformed KZ({block_start}:{block_end}) block")
        actual = values[index - block_start]
        if actual != expected:
            fail(f"{case_name}: KZ({index})={actual}, expected {expected}")

    dat10 = parse_numeric_rows(case_dir / "dat.10", 7)
    datsev_path = case_dir / "datsev.21"
    datsev_lines = [
        line for line in datsev_path.read_text().splitlines() if line.strip()
    ]
    if not datsev_lines:
        fail(f"{case_name}: datsev.21 is empty")
    if len(datsev_lines[0].split()) != 5:
        fail(f"{case_name}: datsev.21 header must have five columns")
    try:
        header = [float(field) for field in datsev_lines[0].split()]
    except ValueError as exc:
        fail(f"{case_name}: datsev.21 header is not numeric: {exc}")
    if not all(math.isfinite(value) for value in header):
        fail(f"{case_name}: datsev.21 header contains a non-finite value")
    stars: list[list[float]] = []
    for line_number, raw in enumerate(datsev_lines[1:], 2):
        fields = raw.split()
        if len(fields) != 6:
            fail(
                f"{datsev_path}: line {line_number} has {len(fields)} columns, expected 6"
            )
        try:
            row = [float(field) for field in fields]
        except ValueError as exc:
            fail(f"{datsev_path}: line {line_number} is not numeric: {exc}")
        if not all(math.isfinite(value) for value in row):
            fail(f"{datsev_path}: line {line_number} contains a non-finite value")
        stars.append(row)

    if len(dat10) != 50:
        fail(f"{case_name}: dat.10 has {len(dat10)} bodies, expected 50")
    if len(stars) != 50:
        fail(f"{case_name}: datsev.21 has {len(stars)} stars, expected 50")

    if int(header[1]) != 50 or header[1] != 50:
        fail(f"{case_name}: datsev.21 NZERO={header[1]}, expected 50")

    mean_mass = header[3]
    total_mass = 50.0 * mean_mass
    normalized_mass = sum(row[0] for row in dat10)
    if not math.isclose(normalized_mass, 1.0, rel_tol=0.0, abs_tol=2e-13):
        fail(f"{case_name}: normalized BODY mass sums to {normalized_mass}, expected 1")

    zmb = namelist_number(inp, "ZMBAR")
    rbar = namelist_number(inp, "RBAR")
    if not math.isclose(zmb, mean_mass, rel_tol=0.0, abs_tol=2e-13):
        fail(f"{case_name}: ZMBAR={zmb} but datsev mean mass={mean_mass}")
    if not math.isclose(rbar, header[2], rel_tol=0.0, abs_tol=2e-15):
        fail(f"{case_name}: RBAR={rbar} but datsev RBAR={header[2]}")

    for body_index, (body, star) in enumerate(zip(dat10, stars), 1):
        mass_from_body = body[0] * total_mass
        current_mass = star[0]
        if not math.isclose(
            mass_from_body, current_mass, rel_tol=2e-13, abs_tol=2e-13
        ):
            fail(
                f"{case_name}: body {body_index} mass is {mass_from_body} Msun "
                f"from dat.10 but {current_mass} Msun in datsev.21"
            )
        if star[1] != int(star[1]):
            fail(f"{case_name}: body {body_index} has non-integral KSTAR={star[1]}")

    print(f"{case_name}: PASS")


def main() -> None:
    manifest = json.loads((ROOT / "manifest.json").read_text())
    executable_hash = manifest.get("validated_executable_sha256")
    if not isinstance(executable_hash, str) or not re.fullmatch(
        r"[0-9a-f]{64}", executable_hash
    ):
        fail("validated_executable_sha256 must be 64 lowercase hex characters")
    fixtures = manifest.get("fixtures")
    if not isinstance(fixtures, dict) or len(fixtures) != 5:
        fail("manifest must define exactly five fixtures")
    for case_name, case in fixtures.items():
        if not isinstance(case, dict):
            fail(f"{case_name}: malformed manifest entry")
        verify_case(case_name, case)
    print("All five fixture input sets are intact; no simulation was launched.")


if __name__ == "__main__":
    main()
