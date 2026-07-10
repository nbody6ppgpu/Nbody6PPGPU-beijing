#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

make_fake() {
  prefix=$1 compiler=$2
  mkdir -p "$prefix/bin" "$prefix/lib" "$prefix/lib64"
  cp "$tmp/fake-compiler" "$prefix/bin/$compiler"
  chmod +x "$prefix/bin/$compiler"
}

cat > "$tmp/fake-compiler" <<'EOF'
#!/bin/sh
out=
while test $# -gt 0; do
  if test "x$1" = x-o; then shift; out=$1; fi
  shift
done
test -z "$out" || : > "$out"
exit 0
EOF

make_fake "$tmp/cuda" nvcc
make_fake "$tmp/hip" hipcc

configure_case() {
  name=$1 expected=$2
  shift 2
  rm -f "$root/config.status" "$root/config.log" "$root/build/Makefile"
  (cd "$root" && PATH="$tmp/base:$PATH" "$@" \
    --build=x86_64-pc-linux-gnu --host=x86_64-pc-linux-gnu \
    --enable-mcmodel=no > "$tmp/$name.log" 2>&1)
  grep "GPU backend: *$expected" "$tmp/$name.log" >/dev/null
}

mkdir -p "$tmp/base"
configure_case disabled none ./configure --disable-gpu --disable-mpi --enable-simd=no
configure_case cuda cuda env NVCC="$tmp/cuda/bin/nvcc" ./configure --disable-mpi --enable-simd=no
grep '^EXTRAOBJ=.*CUDA_OBJECTS' "$root/build/Makefile" >/dev/null
! grep '^EXTRAOBJ=.*HIP_OBJECTS' "$root/build/Makefile" >/dev/null
configure_case hip hip env HIPCC="$tmp/hip/bin/hipcc" ./configure --with-gpu-backend=hip --disable-mpi --enable-simd=no
grep '^EXTRAOBJ=.*HIP_OBJECTS' "$root/build/Makefile" >/dev/null
! grep '^EXTRAOBJ=.*CUDA_OBJECTS' "$root/build/Makefile" >/dev/null
configure_case priority cuda env NVCC="$tmp/cuda/bin/nvcc" HIPCC="$tmp/hip/bin/hipcc" ./configure --disable-mpi --enable-simd=no
configure_case cuda_prefix cuda ./configure --with-gpu-backend=cuda --with-cuda="$tmp/cuda" --disable-mpi --enable-simd=no
configure_case hip_prefix hip ./configure --with-gpu-backend=hip --with-hip="$tmp/hip" --disable-mpi --enable-simd=no

if (cd "$root" && ./configure --build=x86_64-pc-linux-gnu --host=x86_64-pc-linux-gnu --enable-mcmodel=no --disable-gpu --with-gpu-backend=hip >/dev/null 2>&1); then exit 1; fi
if (cd "$root" && ./configure --build=x86_64-pc-linux-gnu --host=x86_64-pc-linux-gnu --enable-mcmodel=no --with-gpu-backend=bogus >/dev/null 2>&1); then exit 1; fi
if (cd "$root" && PATH="$tmp/base" ./configure --build=x86_64-pc-linux-gnu --host=x86_64-pc-linux-gnu --enable-mcmodel=no --disable-mpi --enable-simd=no >/dev/null 2>&1); then exit 1; fi

echo "GPU configure matrix passed"
