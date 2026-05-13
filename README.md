This is a version for running on ARM-based computers, which has been tested on Oracle Cloud Instance (VPS) with 3.0 GHz Ampere® Altra™ CPU (4 × Neoverse N1 cores) and also the Jülich JUPITER Booster supercomputer with [NVIDIA GH200](https://apps.fz-juelich.de/jsc/hps/jupiter/configuration.html#jupiter-hardware-overview) Superchips.

This version is 99% the same as the standard version, with only two additional steps for configuration before installation.

1. get the SIMDe library: `git clone --depth 1 https://github.com/simd-everywhere/simde`
2. during the configuration, specify the path to the SIMDe library: `./configure --with-simde=[PATH_TO_SIMDe_ROOT] [other configure options]` 

All other documentation is the same as the standard version. [Read the README there](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/tree/dev)
