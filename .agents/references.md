# Documentation Resources and References

### Official documentation
- **Manual source**: `doc/` (git submodule pointing to [nbody6ppgpu/nb6-manual](https://github.com/nbody6ppgpu/nb6-manual), pinned commit; this directory is empty until `git submodule update --init --recursive` is run)
- **Manual PDF**: https://nbody6ppgpu.github.io/nb6-manual-pdf/latest.pdf
- **Overleaf (entry point for manual editing by hand, automatically synced both ways with the manual repo)**: [Overleaf Manual](https://www.overleaf.com/read/hcmxcyffjkzq)
- **README**: [GitHub README](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/blob/stable/README.md)
- **Discussions**: [GitHub Discussions](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/discussions)

### Editing the manual

If you change user-visible behaviour, parameters or IO formats, first search `doc/` with `rg` to check whether the manual needs a matching update. If `doc/` is an empty directory,
first run `git submodule update --init --recursive doc` to fetch the manual content; do not conclude "there is no manual to edit" just because the directory is empty.

The manual source itself belongs to the separate `nb6-manual` repo (automatically synced both ways with Overleaf); do not edit file contents under `doc/` directly in this repo;
after editing, the change is complete only once you have run `cd doc && git add -A && git commit && git push` to the `nb6-manual` origin -- merely
committing the gitlink pointer in the parent repo does **not** persist source changes in `doc/`; the commit inside `doc/` must be pushed to the remote first.
After the push succeeds, return to the root of this repo and run `git add doc && git commit` to advance the pointer to the new commit.

### References

| Topic | Reference | NASA ADS |
|------|------|----------|
| NBODY1-6 review | Aarseth 1999 | [1999PASP..111.1333A](https://ui.adsabs.harvard.edu/abs/1999PASP..111.1333A/) |
| NBODY6++ | Spurzem 1999 | [1999JCoAM.109..407S](https://ui.adsabs.harvard.edu/abs/1999JCoAM.109..407S/) |
| SSE/BSE | Hurley et al. 2005 | [2005MNRAS.363..293H](https://ui.adsabs.harvard.edu/abs/2005MNRAS.363..293H/) |
| NBODY6GPU | Nitadori & Aarseth 2012 | [2012MNRAS.424..545N](https://ui.adsabs.harvard.edu/abs/2012MNRAS.424..545N/) |
| NBODY6++GPU | Wang et al. 2015 | [2015MNRAS.450.4070W](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.4070W/) |
| Stellar evolution update | Kamlah et al. 2022 | [2022MNRAS.511.4060K](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4060K/) |

### Related tools
- **dragon3_pipeline**: [Complete, extensive Python scripts for data analysis (under construction)](https://github.com/kaiwu-astro/dragon3_pipeline)
- **McLuster**: [Initial conditions generator](https://github.com/agostinolev/mcluster)
- **Galpy**: [Milky Way potential library](https://github.com/jobovy/galpy)
