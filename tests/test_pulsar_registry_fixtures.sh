#!/bin/sh
# Static validation of the pulsar-registry end-to-end fixtures in
# tests/pulsar_registry_fixtures/ (checksums, input counts, finite
# values, normalized BODY masses, and dat.10/datsev.21/namelist
# agreement). This does not build or run NBODY6++GPU; it only checks
# that the committed fixture inputs are intact, so it runs in seconds
# and is suitable for CI. This script does not build or execute
# NBODY6++GPU and therefore does not check any of manifest.json's
# `expected` simulation results (final NS counts, ERRTOT, event paths,
# etc.) -- see tests/pulsar_registry_fixtures/README.md, including its
# Reproducibility section, for how to actually run a fixture and for
# known environment sensitivity of the provisional_ns_cleanup case.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

cd "$root"

fixtures_dir="tests/pulsar_registry_fixtures"
test -d "$fixtures_dir" || {
  echo "FAIL: $fixtures_dir not found"
  exit 1
}

if ! python3 "$fixtures_dir/verify_inputs.py" > "$tmp/verify.log" 2>&1; then
  echo "FAIL: verify_inputs.py reported a problem with the fixture inputs"
  cat "$tmp/verify.log"
  exit 1
fi

cat "$tmp/verify.log"

if ! grep -q "All five fixture input sets are intact" "$tmp/verify.log"; then
  echo "FAIL: verify_inputs.py did not report all five fixture sets intact"
  exit 1
fi

echo "Pulsar registry fixture inputs verified."
