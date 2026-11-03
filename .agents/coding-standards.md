# Coding Standards

### Fortran coding standards

- **Format**: Fortran 77 fixed form (6-column indent), or free-form Fortran 90 in a very few files
- **Case**: Mixed (legacy code is upper case, new code may use lower case)
- **Comments**: Start with `C` or `!`
- **Common blocks**: Variables are shared through `COMMON` blocks (see `include/common6.h`)
- **Preprocessor**: Use the `.F` extension to enable preprocessor directives (e.g. `#ifdef`)

### Naming conventions

| Type | Convention | Example |
|------|------|------|
| Subroutines | Lower case, descriptive names | `intgrt`, `ksinit`, `output` |
| Functions | Upper case, short | `ENERGY`, `LAGR` |
| Variables | Upper case (common blocks), lower case (local) | `TIME`, `TPHYS`, `body(i)` |
| Parameters | Upper case | `NMAX`, `KMAX` |
