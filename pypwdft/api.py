# -*- coding: utf-8 -*-

# 
# This file is part of the PyPWDFT distribution 
# Copyright (c) 2024 Ivo Filot
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from importlib import resources
from pathlib import Path

import numpy as np
from mendeleev import element

from .gth import GTHPseudopotential
from .psystem import PeriodicSystem
from .pypwdft import PyPWDFT


_ANGSTROM_TO_BOHR = 1.8897259886
_GTH_FAMILY_FOR_XC = {"lda": "pade", "pbe": "pbe"}


class Structure:
    """
    Atomic structure in a cubic simulation cell.

    Parameters
    ----------
    symbols : sequence of str
        Element symbols in the same order as ``positions``.
    positions : array_like, shape (natom, 3)
        Cartesian atomic coordinates.
    cell : float
        Cubic cell edge length in the selected units.
    units : {'bohr', 'angstrom'}, optional
        Units used by ``positions`` and ``cell``. Internally, coordinates are
        stored in bohr.
    center : bool, optional
        Centre the structure's bounding box in the simulation cell.
    """

    def __init__(self, symbols, positions, cell, *, units="bohr", center=True):
        """
        Initialize a Structure object with atomic symbols, positions, and cell
        size.
        """
        # build cubic cell and validate inputs
        if units not in {"bohr", "angstrom"}:
            raise ValueError("units must be 'bohr' or 'angstrom'.")
        scale = 1.0 if units == "bohr" else _ANGSTROM_TO_BOHR
        cell = float(cell) * scale
        if not np.isfinite(cell) or cell <= 0:
            raise ValueError("cell must be positive and finite.")

        # validate symbols and positions
        normalised_symbols = tuple(element(symbol).symbol for symbol in symbols)
        coordinates = np.asarray(positions, dtype=float)
        if coordinates.shape != (len(normalised_symbols), 3):
            raise ValueError("positions must have shape (len(symbols), 3).")
        if len(coordinates) == 0:
            raise ValueError("a structure must contain at least one atom.")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("positions must contain only finite values.")
        coordinates = coordinates * scale

        # optionally center the structure in the simulation cell
        if center:
            midpoint = (coordinates.min(axis=0) + coordinates.max(axis=0)) / 2.0
            coordinates = coordinates + cell / 2.0 - midpoint
        if np.any(coordinates < 0) or np.any(coordinates >= cell):
            raise ValueError(
                "all atoms must lie inside the cell; use center=True or a "
                "larger cell."
            )

        self._symbols = normalised_symbols
        self._positions = coordinates.copy()
        self._cell = cell

    @property
    def symbols(self):
        """Element symbols as an immutable tuple."""
        return self._symbols

    @property
    def positions(self):
        """A copy of the Cartesian positions in bohr."""
        return self._positions.copy()

    @property
    def cell(self):
        """Cubic cell edge length in bohr."""
        return self._cell

    @property
    def atomic_numbers(self):
        """Atomic numbers in structure order."""
        return np.asarray([element(symbol).atomic_number for symbol in self.symbols])

    @classmethod
    def from_name(cls, name, *, cell=10.0, center=True):
        """Load one of PyPWDFT's bundled XYZ structures."""
        filename = resources.files("pypwdft").joinpath(
            "molecules", f"{str(name).lower()}.xyz"
        )
        if not filename.is_file():
            raise ValueError(f"unknown bundled structure: {name!r}")
        with resources.as_file(filename) as path:
            return cls.from_xyz(path, cell=cell, center=center)

    @classmethod
    def from_xyz(cls, path, *, cell=10.0, center=True):
        """Load an XYZ file whose coordinates are expressed in angstrom."""
        path = Path(path)
        with path.open(encoding="utf-8") as handle:
            lines = handle.readlines()
        try:
            natoms = int(lines[0].strip())
        except (IndexError, ValueError) as exc:
            raise ValueError(f"invalid XYZ file: {path}") from exc
        atom_lines = lines[2:2 + natoms]
        if len(atom_lines) != natoms:
            raise ValueError(f"invalid XYZ file: expected {natoms} atoms")

        symbols = []
        positions = []
        for line in atom_lines:
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"invalid XYZ atom line: {line.rstrip()!r}")
            symbols.append(fields[0])
            try:
                positions.append([float(value) for value in fields[1:4]])
            except ValueError as exc:
                raise ValueError(f"invalid XYZ atom line: {line.rstrip()!r}") from exc
        positions = np.asarray(positions) * _ANGSTROM_TO_BOHR
        return cls(symbols, positions, cell, units="bohr", center=center)

    def _to_periodic_system(self, cutoff, *, density_cutoff=None):
        """Build the internal numerical system used by the solver."""
        system = PeriodicSystem(
            self.cell, ecut=cutoff, density_ecut=density_cutoff
        )
        for position, atomic_number in zip(self._positions, self.atomic_numbers):
            system.add_atom(*position, int(atomic_number))
        return system

    def __len__(self):
        return len(self.symbols)

    def __repr__(self):
        formula = "".join(self.symbols)
        return f"Structure({formula}, atoms={len(self)}, cell={self.cell:.3f} bohr)"

#-------------------------------------------------------------------------------

@dataclass(frozen=True)
class GTH:
    """
    Configuration for a GTH pseudopotential.

    ``family=None`` automatically chooses PADE for LDA and PBE for PBE.
    """

    family: str | None = None
    path: str | Path | None = None
    charges: Mapping | None = None

    def _build(self, system, xc):
        family = self.family or _GTH_FAMILY_FOR_XC[xc]
        _validate_gth_family(family, xc, custom_path=self.path is not None)
        return GTHPseudopotential(
            system,
            path=self.path,
            charge_overrides=self.charges,
            family=family,
        )

#-------------------------------------------------------------------------------

@dataclass(frozen=True)
class SCFSettings:
    """Convergence, eigensolver, and density-mixing settings.

    Parameters
    ----------
    convergence : float, optional
        Total-energy convergence threshold in Hartree.
    density_convergence : float, optional
        RMS density-residual threshold. Defaults to ``convergence``.
    max_iterations : int, optional
        Maximum number of SCF iterations.
    bands : int, optional
        Number of orbitals to calculate. Occupied orbitals are always included.
    verbosity : int, optional
        Set to a non-zero value to print the calculation and iteration tables.
    mixing : {'pulay', 'linear'}, optional
        Density-mixing algorithm. Pulay is the default.
    mixing_fraction : float, optional
        Damping applied to each density residual. Must be in ``(0, 1]``.
    mixing_history : int, optional
        Number of recent residuals retained by Pulay mixing.
    """

    convergence: float = 1e-5
    density_convergence: float | None = None
    max_iterations: int = 100
    bands: int | None = None
    verbosity: int = 0
    mixing: str = "pulay"
    mixing_fraction: float = 0.5
    mixing_history: int = 6

    def __post_init__(self):
        if self.convergence <= 0:
            raise ValueError("convergence must be positive.")
        if self.density_convergence is not None and self.density_convergence <= 0:
            raise ValueError("density_convergence must be positive.")
        if self.max_iterations < 1:
            raise ValueError("max_iterations must be at least one.")
        if self.bands is not None and self.bands < 1:
            raise ValueError("bands must be positive.")
        if not isinstance(self.verbosity, (bool, int, np.integer)):
            raise ValueError("verbosity must be an integer.")
        mixing = str(self.mixing).lower()
        if mixing not in {"linear", "pulay"}:
            raise ValueError("mixing must be 'linear' or 'pulay'.")
        object.__setattr__(self, "mixing", mixing)
        if not 0 < self.mixing_fraction <= 1:
            raise ValueError("mixing_fraction must be in the interval (0, 1].")
        if (
            not isinstance(self.mixing_history, (int, np.integer))
            or self.mixing_history < 2
        ):
            raise ValueError("mixing_history must be an integer of at least 2.")

#-------------------------------------------------------------------------------

@dataclass(frozen=True)
class EnergyComponents:
    """Named energy components in Hartree."""

    total: float
    kinetic: float
    electron_ion: float
    nonlocal_: float
    hartree: float
    xc: float
    ion_ion: float

    @property
    def ionic(self):
        """Combined local, non-local, and ion-ion energy."""
        return self.electron_ion + self.nonlocal_ + self.ion_ion

#-------------------------------------------------------------------------------

@dataclass(frozen=True)
class OrbitalSet:
    """Orbital energies and wavefunctions returned by an SCF calculation."""

    energies: np.ndarray
    real_space: np.ndarray
    reciprocal_space: np.ndarray

#-------------------------------------------------------------------------------

@dataclass(frozen=True)
class BasisInfo:
    """Plane-wave basis and real-space grid information."""

    cutoff: float
    density_cutoff: float
    plane_waves: int
    wavefunction_grid: int
    density_grid: int

#-------------------------------------------------------------------------------

@dataclass(frozen=True)
class SCFInfo:
    """SCF convergence and timing information."""

    converged: bool
    iterations: int
    energy_residual: float
    density_residual: float
    elapsed_time: float
    backend: str
    mixing: str
    mixing_fraction: float
    mixing_history: int

#-------------------------------------------------------------------------------

class DFTResult:
    """Typed result of a completed DFT calculation."""

    def __init__(self, data, *, system, noccupied):
        self._system = system
        self._noccupied = int(noccupied)
        self._density = np.asarray(data["edens"])
        self.energy = EnergyComponents(
            total=float(data["Etot"]),
            kinetic=float(data["Ekin"]),
            electron_ion=float(data["Enuc"]),
            nonlocal_=float(data.get("Enonloc", 0.0)),
            hartree=float(data["Erep"]),
            xc=float(data["Exc"]),
            ion_ion=float(data["Eewald"]),
        )
        self.orbitals = OrbitalSet(
            energies=np.asarray(data["orbe"]),
            real_space=np.asarray(data["orbc_rs"]),
            reciprocal_space=np.asarray(data["orbc_fft"]),
        )
        self.basis = BasisInfo(
            cutoff=float(data["ecut"]),
            density_cutoff=float(data["density_ecut"]),
            plane_waves=int(data["npw"]),
            wavefunction_grid=int(data["wavefunction_npts"]),
            density_grid=int(data["density_npts"]),
        )
        self.scf = SCFInfo(
            converged=bool(data["converged"]),
            iterations=int(data["iterations"]),
            energy_residual=float(data["energy_residual"]),
            density_residual=float(data["density_residual"]),
            elapsed_time=float(data["ttime"]),
            backend=str(data["fft"]),
            mixing=str(data["mixing"]),
            mixing_fraction=float(data["mixing_fraction"]),
            mixing_history=int(data["mixing_history"]),
        )

    @property
    def converged(self):
        return self.scf.converged

    @property
    def iterations(self):
        return self.scf.iterations

    @property
    def density(self):
        return self._density

    @property
    def noccupied(self):
        return self._noccupied

    def plot_orbitals(
        self,
        *,
        occupied=True,
        indices=None,
        plane="xy",
        columns=3,
        save=None,
        show_imaginary_norm=True,
        **plot_options,
    ):
        """Plot selected orbitals and return the Matplotlib figure and axes."""
        from .visualization import ContourPlotter

        if indices is None:
            stop = self.noccupied if occupied else len(self.orbitals.energies)
            indices = np.arange(stop)
        else:
            indices = np.asarray(indices, dtype=int)
        if indices.ndim != 1 or len(indices) == 0:
            raise ValueError("indices must select at least one orbital.")
        if np.any(indices < 0) or np.any(indices >= len(self.orbitals.energies)):
            raise IndexError("orbital index out of range.")
        if not isinstance(columns, (int, np.integer)) or columns < 1:
            raise ValueError("columns must be a positive integer.")

        columns = min(int(columns), len(indices))
        rows = int(np.ceil(len(indices) / columns))
        selected = {
            "orbc_rs": self.orbitals.real_space[indices],
            "orbe": self.orbitals.energies[indices],
        }
        labels = plot_options.pop("labels", None)
        if labels is None:
            labels = [rf"$\psi_{{{index + 1}}}$" for index in indices]
        return ContourPlotter.build_contourplot(
            selected,
            self._system,
            filename=save,
            plane=plane,
            nrows=rows,
            ncols=columns,
            labels=labels,
            show_imaginary_norm=show_imaginary_norm,
            **plot_options,
        )

    def plot_density(self, *, plane="xy", save=None, **plot_options):
        """Plot a central plane through the converged electron density."""
        from .visualization import ContourPlotter

        return ContourPlotter.build_contourplot(
            {"orbc_rs": self.density[np.newaxis, ...]},
            self._system,
            filename=save,
            plane=plane,
            nrows=1,
            ncols=1,
            labels=[r"$\rho$"],
            plot_energies=False,
            **plot_options,
        )

    def __repr__(self):
        return (
            f"DFTResult(total_energy={self.energy.total:.8f} Ha, "
            f"converged={self.converged}, iterations={self.iterations})"
        )

#-------------------------------------------------------------------------------

class PWDFT:
    """
    High-level plane-wave DFT calculation.

    Parameters
    ----------
    structure : Structure
        Atomic structure and cubic simulation cell.
    cutoff : float, optional
        Wavefunction cutoff in Hartree. Required for ``Structure`` input.
    density_cutoff : float, optional
        Density cutoff in Hartree. Defaults to four times ``cutoff``.
    xc : {'lda', 'svwn5', 'pbe'}, optional
        Exchange-correlation functional.
    pseudopotential : {'gth', 'all-electron'}, GTH, object, or None
        Ionic model. The default GTH family automatically matches ``xc``.
    device : {'cpu', 'cuda'}, optional
        Compute device. CPU is the predictable default.
    fft_backend : {'pyfftw', 'numpy', 'scipy', 'cupy'}, optional
        Expert override for the numerical backend.
    """

    def __init__(
        self,
        structure,
        *,
        cutoff=None,
        density_cutoff=None,
        xc="pbe",
        pseudopotential="gth",
        device="cpu",
        fft_backend=None,
    ):
        # validate inputs and build the internal PyPWDFT engine
        xc = str(xc).lower()
        if xc not in _GTH_FAMILY_FOR_XC:
            raise ValueError("xc must be 'lda' or 'pbe'.")

        # build the internal numerical system used by the solver
        if not isinstance(structure, Structure):
            raise TypeError("structure must be a Structure.")
        if cutoff is None:
            raise ValueError("cutoff is required.")
        system = structure._to_periodic_system(
            cutoff, density_cutoff=density_cutoff
        )

        # validate device and FFT backend
        device = str(device).lower()
        if device not in {"cpu", "cuda"}:
            raise ValueError("device must be 'cpu' or 'cuda'.")
        if fft_backend is None:
            fft_backend = "cupy" if device == "cuda" else "pyfftw"
        fft_backend = str(fft_backend).lower()
        if device == "cuda" and fft_backend != "cupy":
            raise ValueError("device='cuda' requires fft_backend='cupy'.")
        if device == "cpu" and fft_backend == "cupy":
            raise ValueError("fft_backend='cupy' requires device='cuda'.")

        # build the ionic model and store all user-facing attributes
        ionic_model = _build_pseudopotential(
            pseudopotential, system=system, xc=xc
        )
        self.structure = structure
        self._system = system
        self.cutoff = system.get_ecut()
        self.density_cutoff = system.get_density_ecut()
        self.xc = xc
        self._ionic_model = ionic_model
        self.ionic_model = (
            "all-electron"
            if ionic_model is None
            else f"gth-{ionic_model.family}"
        )
        self.device = device
        self.fft_backend = fft_backend
        self._engine = PyPWDFT(
            system,
            fft=fft_backend,
            pseudopotential=ionic_model,
            functional=xc,
        )

    def run(
        self,
        settings=None,
        *,
        convergence=None,
        density_convergence=None,
        max_iterations=None,
        bands=None,
        verbosity=None,
        mixing=None,
        mixing_fraction=None,
        mixing_history=None,
    ):
        """Run the SCF calculation and return a :class:`DFTResult`.

        Keyword arguments override the corresponding values in ``settings``.
        In most calculations the default Pulay mixer is preferable to linear
        mixing. If convergence oscillates, reduce ``mixing_fraction``.
        """
        # validate settings and override with any explicit arguments
        if settings is None:
            settings = SCFSettings()
        elif not isinstance(settings, SCFSettings):
            raise TypeError("settings must be an SCFSettings instance.")

        # override settings with any explicit arguments
        overrides = {
            "convergence": convergence,
            "density_convergence": density_convergence,
            "max_iterations": max_iterations,
            "bands": bands,
            "verbosity": verbosity,
            "mixing": mixing,
            "mixing_fraction": mixing_fraction,
            "mixing_history": mixing_history,
        }
        settings = replace(
            settings,
            **{key: value for key, value in overrides.items() if value is not None},
        )
        data = self._engine.scf(
            tol=settings.convergence,
            density_tol=settings.density_convergence,
            maxiter=settings.max_iterations,
            nsol=settings.bands,
            verbose=bool(settings.verbosity),
            mixing=settings.mixing,
            mixing_fraction=settings.mixing_fraction,
            mixing_history=settings.mixing_history,
        )
        nelectrons = (
            self._ionic_model.nelec
            if self._ionic_model is not None
            else self._system.get_nelec()
        )
        return DFTResult(
            data, system=self._system, noccupied=int(nelectrons) // 2
        )

    def __repr__(self):
        ionic = (
            "all-electron"
            if self._ionic_model is None
            else f"GTH-{self._ionic_model.family.upper()}"
        )
        return (
            f"PWDFT(xc={self.xc!r}, cutoff={self.cutoff:g} Ha, "
            f"ionic_model={ionic!r}, device={self.device!r})"
        )


def _validate_gth_family(family, xc, *, custom_path=False):
    """
    Validate that the GTH family matches the selected exchange-correlation
    functional. If a custom path is provided, the family is not validated.
    """
    family = str(family).lower()
    if not custom_path and family != _GTH_FAMILY_FOR_XC[xc]:
        raise ValueError(
            f"GTH family {family!r} does not match xc={xc!r}; use "
            f"{_GTH_FAMILY_FOR_XC[xc]!r}."
        )


def _build_pseudopotential(specification, *, system, xc):
    """
    Build a pseudopotential from the given specification.
    """
    if specification is None:
        return None
    if isinstance(specification, GTH):
        return specification._build(system, xc)
    if isinstance(specification, str):
        name = specification.lower()
        if name == "all-electron":
            return None
        if name == "gth":
            return GTH()._build(system, xc)
        if name in {"gth-pbe", "gth-pade"}:
            return GTH(family=name.removeprefix("gth-"))._build(system, xc)
        raise ValueError(
            "pseudopotential must be 'gth', 'gth-pbe', 'gth-pade', "
            "'all-electron', or a GTH configuration."
        )
    raise TypeError("unsupported pseudopotential specification.")
