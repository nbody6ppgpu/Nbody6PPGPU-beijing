#!/bin/sh
# Regenerate ./configure from configure.ac with the pinned autoconf toolchain.
#
#   tool/regen-configure.sh           regenerate configure in place
#   tool/regen-configure.sh --check   fail if the committed configure is stale
#
# Always use this instead of running autoconf directly: see
# tool/autoconf-env/Dockerfile for why the exact toolchain matters.
#
# --check always requires docker, because its whole job is to reproduce CI's
# verdict.  Plain regeneration accepts --allow-system-autoconf as an escape
# hatch on a machine without docker; the result is then only as good as whatever
# autoconf that machine has, and CI still decides.
set -eu

repo_root=$(CDPATH= cd "$(dirname "$0")/.." && pwd)
image=nbody6ppgpu-autoconf:2.71-2
want_autoconf=2.71-2
want_m4=1.4.18-5ubuntu2
check=no
allow_system=no

for arg in "$@"; do
    case $arg in
        --check) check=yes ;;
        --allow-system-autoconf) allow_system=yes ;;
        *) echo "regen-configure: unknown option: $arg" >&2; exit 2 ;;
    esac
done

if [ "$check" = yes ] && [ "$allow_system" = yes ]; then
    echo "regen-configure: --check cannot be combined with --allow-system-autoconf;" >&2
    echo "regen-configure: a check that does not use the pinned toolchain proves nothing." >&2
    exit 2
fi

assert_pinned_versions() {
    got_autoconf=$(docker run --rm "$image" dpkg-query -W -f='${Version}' autoconf)
    got_m4=$(docker run --rm "$image" dpkg-query -W -f='${Version}' m4)
    if [ "$got_autoconf" != "$want_autoconf" ] || [ "$got_m4" != "$want_m4" ]; then
        echo "regen-configure: pinned image drifted: autoconf $got_autoconf (want $want_autoconf)," >&2
        echo "regen-configure: m4 $got_m4 (want $want_m4).  Refusing to generate." >&2
        exit 1
    fi
}

regen() {
    if command -v docker >/dev/null 2>&1; then
        docker build -q -t "$image" "$repo_root/tool/autoconf-env" >/dev/null
        assert_pinned_versions
        echo "regen-configure: using pinned autoconf $got_autoconf / m4 $got_m4" >&2
        # The source tree is mounted read-only and autoconf runs from a throwaway
        # working directory, so the repo's autom4te.cache is neither read nor
        # written.  That matters: a cache left by a different autoconf silently
        # poisons the output, and autoconf also short-circuits when it thinks the
        # output is current, which would let a hand-edited configure pass --check.
        docker run --rm -u "$(id -u):$(id -g)" \
            -v "$repo_root:/src:ro" -v "$outdir:/out" \
            --tmpfs /work:rw,mode=1777 -w /work "$image" \
            autoconf --force -I /src -o "/out/$outbase" /src/configure.ac
    elif [ "$allow_system" = yes ]; then
        echo "regen-configure: WARNING using the system autoconf, not the pinned one:" >&2
        autoconf --version | sed -n '1s/^/regen-configure:   /p' >&2
        echo "regen-configure: CI compares against autoconf $want_autoconf; expect disagreement." >&2
        ( cd "$repo_root" && autoconf --force -o "$out" )
    else
        cat >&2 <<MSG
regen-configure: docker is required so that everyone regenerates configure with
the same toolchain (autoconf $want_autoconf, m4 $want_m4).  Install docker, or -- for
plain regeneration only -- pass --allow-system-autoconf to use the local one.
MSG
        exit 1
    fi
}

outdir=$(mktemp -d)
trap 'rm -rf "$outdir"' 0
outbase=configure.regen
out=$outdir/$outbase

regen

if [ "$check" = yes ]; then
    if ! diff -u "$repo_root/configure" "$out"; then
        cat >&2 <<'MSG'

configure is out of sync with configure.ac (diff above: - committed, + regenerated).
Fix it by running

    tool/regen-configure.sh

and committing the result.  Note that dnl strips itself to end of line but not
the whitespace before it, so keep dnl comments in configure.ac at column 0.
MSG
        exit 1
    fi
    echo "configure is in sync with configure.ac"
else
    cat "$out" > "$repo_root/configure"
    chmod +x "$repo_root/configure"
    echo "regenerated $repo_root/configure"
fi
