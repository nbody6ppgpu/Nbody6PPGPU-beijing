# Pinned toolchain for regenerating ./configure from configure.ac.
#
# ./configure is committed (users on clusters run it without autoconf installed),
# and CI enforces that it is byte-identical to what configure.ac generates.  That
# only works if everyone regenerates it with the same toolchain: autoconf 2.71-2
# and 2.71-3 differ in the cache variable names they emit (ac_cv_prog_cxx_11 vs
# ac_cv_prog_cxx_cxx11), so "autoconf 2.71" on the developer's distro is not a
# tight enough contract.
#
# The base is pinned by the multi-arch index digest of ubuntu:22.04, so this
# builds identically on x86_64 and arm64 runners.  autoconf and m4 -- the two
# packages whose versions actually shape the generated script -- are pinned by
# Debian revision; tool/regen-configure.sh asserts both at run time, so a silent
# archive change turns into a loud failure rather than a mystery diff.  perl is
# deliberately left floating: autom4te runs on it but does not vary its output
# with it, and pinning it would break the build on every jammy security update.
#
# This is reproducible enough to be trustworthy, but it is not hermetic: the apt
# indexes it reads are mutable.  Making it fully hermetic means publishing this
# image to a registry once and pinning it by its own digest, which needs write
# access to the org's package registry.
FROM ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc

RUN apt-get update -qq \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      autoconf=2.71-2 \
      m4=1.4.18-5ubuntu2 \
 && rm -rf /var/lib/apt/lists/*
