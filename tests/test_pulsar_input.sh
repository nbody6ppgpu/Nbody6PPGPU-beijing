#!/bin/sh
# Accept/reject tests for the KZ(29)/INPULSAR input contract (READPULSAR
# in input.F). Builds a minimal serial, no-GPU binary once and checks
# that each fixture in tests/pulsar_input/ is accepted or rejected as
# expected, per merge-preparation/merge-plan.md section 10.3.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

cd "$root"
rm -f config.status config.log build/Makefile include/params.h
./configure --with-par=1k --disable-gpu --disable-mpi > "$tmp/configure.log" 2>&1
make -C build clean > /dev/null 2>&1 || true
make clean > /dev/null 2>&1 || true
make > "$tmp/build.log" 2>&1
bin="$root/build/nbody6++.avx.hdf5"
test -x "$bin" || bin=$(ls "$root"/build/nbody6++.avx*.hdf5 2>/dev/null | head -1)
test -x "$bin"

export OMP_STACKSIZE=1024M
export OMP_NUM_THREADS=2
ulimit -s unlimited 2>/dev/null || true

run_case() {
  name=$1
  rundir="$tmp/run_$name"
  mkdir -p "$rundir"
  (cd "$rundir" && timeout 120 "$bin" < "$root/tests/pulsar_input/$name.inp" \
    > "$tmp/$name.log" 2>&1) || true
}

fail=0

run_case reject_kz29_out_of_range
if ! grep -q "FATAL ERROR: KZ(29) must be 0, 1, or 2" "$tmp/reject_kz29_out_of_range.log"; then
  echo "FAIL: reject_kz29_out_of_range did not report the expected FATAL ERROR"
  fail=1
fi

run_case reject_missing_inpulsar
if ! grep -q "FATAL ERROR: missing or invalid INPULSAR" "$tmp/reject_missing_inpulsar.log"; then
  echo "FAIL: reject_missing_inpulsar did not report the expected FATAL ERROR"
  fail=1
fi

run_case reject_acc_ce_nonzero
if ! grep -q "FATAL ERROR: PSR_ACC_CE must be 0" "$tmp/reject_acc_ce_nonzero.log"; then
  echo "FAIL: reject_acc_ce_nonzero did not report the expected FATAL ERROR"
  fail=1
fi

run_case accept_kz29_1_valid
if grep -q "FATAL ERROR" "$tmp/accept_kz29_1_valid.log"; then
  echo "FAIL: accept_kz29_1_valid unexpectedly reported a FATAL ERROR"
  fail=1
fi
if ! grep -q "END RUN" "$tmp/accept_kz29_1_valid.log"; then
  echo "FAIL: accept_kz29_1_valid did not reach END RUN"
  fail=1
fi
if ! grep -q "PSR PARAMS" "$tmp/accept_kz29_1_valid.log"; then
  echo "FAIL: accept_kz29_1_valid did not print the parsed PSR PARAMS"
  fail=1
fi

if [ "$fail" -ne 0 ]; then
  echo "One or more pulsar input tests failed; logs kept would be in $tmp (removed on exit)"
  exit 1
fi

echo "All pulsar input accept/reject tests passed."
