#!/bin/sh
# Standalone invariant checks for PULSARINIT/PULSAREVO (pulsar.F),
# per merge-preparation/merge-plan.md section 10.4. No lifecycle hook
# calls these routines yet, so this compiles a small driver
# (tests/pulsar_physics/test_pulsar_invariants.f) directly against
# the built pulsar.o/ran2.o and calls them itself.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

cd "$root"
rm -f config.status config.log build/Makefile include/params.h
./configure --with-par=1k --disable-gpu --disable-mpi > "$tmp/configure.log" 2>&1
make clean > /dev/null 2>&1 || true
make > "$tmp/build.log" 2>&1

fc=$(sed -n 's/^FC[[:space:]]*=[[:space:]]*//p' build/Makefile | head -1)
test -n "$fc" || fc=gfortran

"$fc" -O2 -I extra_inc/nompi -I include -fPIC -mcmodel=large \
  -o "$tmp/test_pulsar_invariants" \
  tests/pulsar_physics/test_pulsar_invariants.f \
  build/pulsar.o build/ran2.o \
  > "$tmp/link.log" 2>&1 || {
    echo "FAIL: could not build test_pulsar_invariants"
    cat "$tmp/link.log"
    exit 1
  }

"$tmp/test_pulsar_invariants"
