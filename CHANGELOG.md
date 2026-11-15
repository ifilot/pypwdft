# Changelog

All notable changes to PyPWDFT are documented in this file. The project
follows [Semantic Versioning](https://semver.org/).

## [0.3.0] - 2026-11-15

This release replaces the public interface. Scripts written for 0.2.x need to
be ported; see *Migrating from 0.2.x* below.

### Added

- High-level API: `Structure`, `PWDFT`, `SCFSettings`, `GTH`, and `DFTResult`
  with typed `EnergyComponents`, `OrbitalSet`, `BasisInfo`, and `SCFInfo`
  results.
- GTH frozen-core pseudopotentials with local and non-local terms. The
  GTH-PADE and GTH-PBE families are bundled (from eminus / CP2K) and matched
  automatically to the chosen functional; `GTH(charges=...)` selects
  alternative valence configurations.
- PBE exchange-correlation functional alongside LDA (Slater + VWN5).
- Spherical plane-wave cutoff and a separate density cutoff
  (`density_cutoff`, default four times `cutoff`).
- Block Davidson eigensolver, restarted from the previous SCF orbitals, as the
  default. The ARPACK Lanczos solver remains available through
  `eigensolver="arpack"` and is checked against missed degenerate states.
- Pulay density mixing (default) and linear mixing; convergence is judged on
  both the energy change and the density residual.
- Optional CuPy/CUDA backend (`device="cuda"`) and a `threads` option that
  avoids BLAS/FFT thread contention.
- Real-valued Γ-point orbitals, with rotation inside degenerate levels.
- Contour plotting through `DFTResult.plot_orbitals` (including
  `plane="auto"`, `ngrid`, `tick_rotation`, `show_imaginary_norm`) and
  `DFTResult.plot_density`.
- Rewritten Sphinx documentation (installation, getting started, background,
  usage, API reference), published on GitHub Pages.
- LaTeX technical note on the pseudopotential implementation
  (`technical-notes/`).
- New examples: H2, CO, benzene, and cutoff-convergence studies.

### Changed

- `matplotlib` and `threadpoolctl` are now runtime dependencies.
- The package version is read from the installed distribution metadata
  (`pypwdft.__version__`); `pypwdft/_version.py` is gone.
- Packaging is defined entirely in `pyproject.toml`; `setup.py` is gone.
- Documentation is built and deployed to GitHub Pages by GitHub Actions.

### Removed

- `PyPWDFT`, `PeriodicSystem`, and `SystemBuilder` are no longer exported from
  `pypwdft`; `SystemBuilder` has been removed entirely.
- The `mendeleev` dependency.
- Conda packages: PyPWDFT is no longer built for the `ifilot` Anaconda
  channel; install it from PyPI with `pip install pypwdft`.

### Migrating from 0.2.x

```python
# 0.2.x
from pypwdft import PyPWDFT, SystemBuilder
system = SystemBuilder().from_name("CH4", npts=32, sz=10)
res = PyPWDFT(system).scf()

# 0.3.0
from pypwdft import PWDFT, Structure
structure = Structure.from_name("CH4", cell=10)
result = PWDFT(structure, cutoff=40, xc="lda", pseudopotential="gth").run()
print(result.energy.total)
```

Calculations are now specified with a kinetic-energy `cutoff` (Hartree)
instead of a number of grid points. 0.2.x was all-electron only; pass
`pseudopotential="all-electron"` to reproduce that model.

## [0.2.2] - 2025-02-11

### Fixed

- Conda packaging: the build workflow now uploads the `.conda` package format
  to the Anaconda channel instead of `.tar.bz2`.

## [0.2.1] - 2025-02-11

No functional changes to the library.

### Added

- Bundled BH3 molecule (`SystemBuilder().from_name("BH3")`).
- Usage documentation on visualizing molecular orbitals: two-dimensional
  projections and isosurfaces, with upsampling by quintic interpolation or in
  the frequency domain.
- Examples for CH4 and CO isosurfaces, orbital orthonormalization, and
  real/reciprocal-space transforms.
- PyPI, Anaconda, and build-status badges on the documentation front page.

### Changed

- CI uses version 4 of the GitHub artifact actions; the documentation
  workflow builds from `master`.

## [0.2.0] - 2024-09-30

No functional changes from 0.1.0.

## [0.1.0] - 2024-09-30

Initial release.

### Added

- Plane-wave DFT solver for isolated molecules in a cubic, periodic unit cell
  with Γ-point sampling (`PyPWDFT`, `PeriodicSystem`).
- All-electron Coulomb potential (no pseudopotentials), with an Ewald sum for
  the nuclear repulsion.
- LDA exchange-correlation: Slater exchange and VWN5 correlation.
- A single basis set and FFT grid for both the orbitals and the electron
  density, specified by the number of grid points per direction.
- Self-consistent field procedure (`PyPWDFT.scf()`) with ARPACK (`eigsh`)
  diagonalization and damped density mixing.
- Choice of FFT module: NumPy, SciPy, or pyFFTW (default).
- `SystemBuilder` with bundled molecules (benzene, BF3, CH4, CO, CO2,
  ethylene, H2, H2O, He, LiH, NH3) and XYZ-file input.
- Sphinx documentation, examples, and wheels on PyPI and the `ifilot`
  Anaconda channel.

[0.3.0]: https://github.com/ifilot/pypwdft/compare/v0.2.2...v0.3.0
[0.2.2]: https://github.com/ifilot/pypwdft/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/ifilot/pypwdft/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/ifilot/pypwdft/releases/tag/v0.2.0
[0.1.0]: https://pypi.org/project/pypwdft/0.1.0/
