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

# GTH formulas and file parsing are adapted from eminus, Copyright 2021 and
# 2023 the eminus developers, licensed under Apache-2.0.

"""Goedecker-Teter-Hutter pseudopotentials for the PyPWDFT basis.

The implementation is intentionally contained in one public class.  The class
owns the pseudopotential parameters, local ionic potential, non-local
projectors, valence-electron count, and ionic energy.  Only an optional
Hamiltonian callback is required from the solver.
"""

import math
from importlib import resources
from pathlib import Path

import numpy as np
from mendeleev import element


class GTHPseudopotential:
    """GTH pseudopotential associated with a :class:`PeriodicSystem`.

    Args:
        system: Periodic system containing positions and atomic numbers.
        path: Directory containing CP2K-style ``Element-qN`` files.  When
            omitted, the bundled LDA/PADE parameter set is used.
        charge_overrides: Optional mapping from element symbols (or atomic
            numbers) to the desired valence charge.  By default the available
            file with the smallest valence charge is selected, matching
            ``eminus``.
        family: Bundled parameter family, ``"pade"`` for LDA or ``"pbe"``
            for PBE. Ignored when ``path`` is supplied.

    Notes:
        The object caches quantities that depend on atomic positions.  Call
        :meth:`rebuild` after modifying the geometry.
    """

    def __init__(self, system, path=None, charge_overrides=None,
                 family="pade"):
        self.system = system
        self.family = family.lower()
        self.path = self._resolve_path(path, self.family)
        self.charge_overrides = charge_overrides or {}
        self.parameters = {}
        self.symbols = self._get_symbols()

        for symbol in dict.fromkeys(self.symbols):
            charge = self.charge_overrides.get(
                symbol,
                self.charge_overrides.get(element(symbol).atomic_number),
            )
            self.parameters[symbol] = self._read_gth(symbol, charge)

        self.valence_charges = np.asarray(
            [self.parameters[symbol]["Zion"] for symbol in self.symbols],
            dtype=int,
        )
        if np.any(self.valence_charges <= 0):
            raise ValueError("Every atom must have a positive GTH valence charge.")

        self.nelec = int(np.sum(self.valence_charges))
        self._local_potential = None
        self._projector_channels = []
        self._device_projector_channels = {}
        self._geometry = None
        self.rebuild()

    @staticmethod
    def _resolve_path(path, family):
        if path is not None:
            result = Path(path)
            if not result.is_dir():
                raise FileNotFoundError(
                    f'GTH pseudopotential directory does not exist: "{result}"'
                )
            return result

        if family not in {"pade", "pbe"}:
            raise ValueError('GTH family must be "pade" or "pbe".')
        package_path = resources.files("pypwdft").joinpath("psp", family)
        if not package_path.is_dir():
            raise FileNotFoundError(
                f"Bundled GTH-{family.upper()} pseudopotentials are not "
                "available; supply path=."
            )
        return package_path

    def _get_symbols(self):
        symbols = []
        for charge in self.system.get_atom_charges():
            atomic_number = int(charge)
            if charge != atomic_number or atomic_number < 1:
                raise ValueError(
                    "GTH atom identities require positive integer atomic numbers."
                )
            symbols.append(element(atomic_number).symbol)
        return tuple(symbols)

    def _find_file(self, symbol, charge):
        if charge is not None:
            candidate = self.path.joinpath(f"{symbol}-q{int(charge)}")
            if not candidate.is_file():
                raise FileNotFoundError(
                    f'No GTH pseudopotential "{candidate.name}" in "{self.path}".'
                )
            return candidate

        candidates = sorted(
            self.path.glob(f"{symbol}-q*"),
            key=lambda item: int(item.name.rsplit("-q", 1)[1]),
        )
        if not candidates:
            raise FileNotFoundError(
                f'No GTH pseudopotential for "{symbol}" in "{self.path}".'
            )
        return candidates[0]

    def _read_gth(self, symbol, charge):
        """Read one CP2K-style GTH parameter file."""
        filename = self._find_file(symbol, charge)
        cloc = np.zeros(4)
        rp = np.zeros(4)
        nproj_l = np.zeros(4, dtype=int)
        h = np.zeros((4, 3, 3))

        with filename.open(encoding="utf-8") as handle:
            handle.readline()  # Element and functional label.
            occupations = handle.readline().split()
            zion = sum(int(value) for value in occupations)

            local = handle.readline().split()
            rloc = float(local[0])
            nlocal = int(local[1])
            values = [float(value) for value in local[2:]]
            if len(values) != nlocal or nlocal > 4:
                raise ValueError(f'Invalid local coefficients in "{filename}".')
            cloc[:nlocal] = values

            lmax = int(handle.readline().split()[0])
            if not 0 <= lmax <= 4:
                raise ValueError(f'Unsupported angular momentum in "{filename}".')

            for l in range(lmax):
                projector = handle.readline().split()
                rp[l] = float(projector[0])
                nproj_l[l] = int(projector[1])
                if not 0 <= nproj_l[l] <= 3:
                    raise ValueError(f'Too many GTH projectors in "{filename}".')

                first_row = [float(value) for value in projector[2:]]
                h[l, 0, :len(first_row)] = first_row
                for row in range(1, nproj_l[l]):
                    values = [float(value) for value in handle.readline().split()]
                    h[l, row, row:row + len(values)] = values

                upper = np.triu(h[l])
                h[l] = upper + np.triu(upper, 1).T

        return {
            "Zion": zion,
            "rloc": rloc,
            "cloc": cloc,
            "lmax": lmax,
            "rp": rp,
            "Nproj_l": nproj_l,
            "h": h,
            "filename": filename,
        }

    @property
    def nprojectors(self):
        """Total number of atom-centred non-local projector functions."""
        return sum(channel[0].shape[1] for channel in self._projector_channels)

    def __getitem__(self, symbol):
        return self.parameters[symbol]

    def rebuild(self):
        """Rebuild all position-dependent local and non-local quantities."""
        positions = np.array(self.system.get_atom_positions(), copy=True)
        if len(positions) != len(self.symbols):
            raise ValueError("Atoms cannot be added after constructing a pseudopotential.")
        self._geometry = positions
        self._local_potential = self._build_local_potential()
        self._projector_channels = self._build_projectors()
        self._device_projector_channels.clear()
        return self

    def _check_geometry(self):
        positions = self.system.get_atom_positions()
        if positions.shape != self._geometry.shape or not np.array_equal(
            positions, self._geometry
        ):
            raise RuntimeError(
                "Atomic positions changed after GTH initialization; call rebuild()."
            )

    def _build_local_potential(self):
        gvec = self.system.get_pw_k().reshape(-1, 3)
        g2 = self.system.get_pw_k2().ravel()
        omega = self.system.get_omega()
        reciprocal_potential = np.zeros(len(g2), dtype=complex)

        for symbol in dict.fromkeys(self.symbols):
            psp = self.parameters[symbol]
            rloc = psp["rloc"]
            c1, c2, c3, c4 = psp["cloc"]
            rloc_g2 = rloc**2 * g2
            exponential = np.exp(-0.5 * rloc_g2)

            with np.errstate(divide="ignore", invalid="ignore"):
                species_potential = (
                    -4 * math.pi * psp["Zion"] * exponential / g2
                    + (2 * math.pi) ** 1.5
                    * rloc**3
                    * exponential
                    * (
                        c1
                        + c2 * (3 - rloc_g2)
                        + c3 * (15 - 10 * rloc_g2 + rloc_g2**2)
                        + c4
                        * (
                            105
                            - 105 * rloc_g2
                            + 21 * rloc_g2**2
                            - rloc_g2**3
                        )
                    )
                )

            zero = g2 == 0
            species_potential[zero] = (
                2 * math.pi * psp["Zion"] * rloc**2
                + (2 * math.pi) ** 1.5
                * rloc**3
                * (c1 + 3 * c2 + 15 * c3 + 105 * c4)
            )
            indices = [i for i, atom_symbol in enumerate(self.symbols)
                       if atom_symbol == symbol]
            structure_factor = np.sum(
                np.exp(1j * gvec @ self._geometry[indices].T), axis=1
            )
            reciprocal_potential += species_potential * structure_factor

        shape = self.system.get_pw_k2().shape
        potential = np.fft.fftn(reciprocal_potential.reshape(shape)) / omega
        return np.real_if_close(potential, tol=1000).real

    def local_potential(self):
        """Return the local ionic potential on the real-space FFT grid."""
        self._check_geometry()
        return self._local_potential

    def _build_projectors(self):
        mask = self.system.get_pw_mask().ravel()
        gvec = self.system.get_pw_k().reshape(-1, 3)[mask]
        gnorm = np.linalg.norm(gvec, axis=1)
        omega = self.system.get_omega()
        channels = []

        for atom_index, symbol in enumerate(self.symbols):
            psp = self.parameters[symbol]
            phase = np.exp(-1j * (gvec @ self._geometry[atom_index]))
            for l in range(psp["lmax"]):
                nproj = int(psp["Nproj_l"][l])
                coupling = np.array(psp["h"][l, :nproj, :nproj], copy=True)
                for m in range(-l, l + 1):
                    beta = np.empty((len(gvec), nproj), dtype=complex)
                    harmonic = self._real_spherical_harmonic(l, m, gvec)
                    for projector in range(nproj):
                        beta[:, projector] = (
                            (-1j) ** l
                            * harmonic
                            * self._evaluate_projector(
                                psp, l, projector + 1, gnorm, omega
                            )
                            * phase
                        )
                    channels.append((beta, coupling))
        return channels

    def apply_nonlocal(self, coefficients):
        """Apply the non-local GTH operator to active PW coefficients."""
        self._check_geometry()
        coefficients = np.asarray(coefficients)
        one_dimensional = coefficients.ndim == 1
        if one_dimensional:
            coefficients = coefficients[:, None]
        if coefficients.ndim != 2 or coefficients.shape[0] != self.system.get_n_plane_waves():
            raise ValueError("Expected active coefficients with shape (npw,) or (npw, nstate).")

        result = np.zeros_like(coefficients, dtype=complex)
        for beta, coupling in self._projector_channels:
            result += beta @ (coupling @ (beta.conj().T @ coefficients))
        return result[:, 0] if one_dimensional else result

    def nonlocal_operator(self, array_module=np):
        """Return a non-local operator using NumPy or a GPU array module.

        GPU projector arrays are copied once and cached until :meth:`rebuild`
        is called.  This keeps projector applications inside the device during
        iterative eigensolver calls.
        """
        self._check_geometry()
        if array_module is np:
            return self.apply_nonlocal

        module_name = array_module.__name__
        try:
            device_id = int(array_module.cuda.runtime.getDevice())
        except AttributeError:
            device_id = None
        cache_key = (module_name, device_id)
        if cache_key not in self._device_projector_channels:
            self._device_projector_channels[cache_key] = [
                (array_module.asarray(beta), array_module.asarray(coupling))
                for beta, coupling in self._projector_channels
            ]
        channels = self._device_projector_channels[cache_key]
        npw = self.system.get_n_plane_waves()

        def apply(coefficients):
            one_dimensional = coefficients.ndim == 1
            if one_dimensional:
                coefficients = coefficients[:, None]
            if coefficients.ndim != 2 or coefficients.shape[0] != npw:
                raise ValueError(
                    "Expected active coefficients with shape (npw,) or "
                    "(npw, nstate)."
                )
            result = array_module.zeros_like(coefficients, dtype=complex)
            for beta, coupling in channels:
                result += beta @ (coupling @ (beta.conj().T @ coefficients))
            return result[:, 0] if one_dimensional else result

        return apply

    def nonlocal_energy(self, occupied_coefficients, occupation=2.0):
        """Return the non-local energy for normalized occupied orbitals."""
        coefficients = np.asarray(occupied_coefficients)
        if coefficients.ndim == 1:
            coefficients = coefficients[:, None]
        applied = self.apply_nonlocal(coefficients)
        return float(
            occupation * np.real(np.sum(coefficients.conj() * applied))
        )

    def ionic_energy(self, gcut=2, gamma=1e-8):
        """Return the periodic ion-ion energy using GTH valence charges."""
        self._check_geometry()
        return self.system.calculate_ewald_sum(
            gcut=gcut, gamma=gamma, charges=self.valence_charges
        )

    @staticmethod
    def _evaluate_projector(psp, l, projector, gnorm, omega):
        radius = psp["rp"][l]
        gr2 = (gnorm * radius) ** 2
        prefactor = 4 * math.pi ** 1.25 * math.sqrt(
            2 ** (l + 1) * radius ** (2 * l + 3) / omega
        )
        value = prefactor * np.exp(-0.5 * gr2)

        if l == 0:
            if projector == 1:
                return value
            if projector == 2:
                return 2 / math.sqrt(15) * (3 - gr2) * value
            if projector == 3:
                return 4 / (3 * math.sqrt(105)) * (15 - 10 * gr2 + gr2**2) * value
        elif l == 1:
            if projector == 1:
                return gnorm / math.sqrt(3) * value
            if projector == 2:
                return 2 * gnorm / math.sqrt(105) * (5 - gr2) * value
            if projector == 3:
                return 4 * gnorm / (3 * math.sqrt(1155)) * (35 - 14 * gr2 + gr2**2) * value
        elif l == 2:
            if projector == 1:
                return gnorm**2 / math.sqrt(15) * value
            if projector == 2:
                return 2 * gnorm**2 / (3 * math.sqrt(105)) * (7 - gr2) * value
        elif l == 3 and projector == 1:
            return gnorm**3 / math.sqrt(105) * value
        raise ValueError(f"No GTH projector for l={l}, projector={projector}.")

    @staticmethod
    def _real_spherical_harmonic(l, m, gvec):
        if l == 0:
            return np.full(len(gvec), 0.5 * math.sqrt(1 / math.pi))

        gnorm = np.linalg.norm(gvec, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            cos_theta = gvec[:, 2] / gnorm
        cos_theta[gnorm < 1e-9] = 0
        sin_theta = np.sqrt(np.maximum(0, 1 - cos_theta**2))
        phi = np.arctan2(gvec[:, 1], gvec[:, 0])

        if l == 1:
            prefactor = 0.5 * math.sqrt(3 / math.pi)
            return {
                -1: prefactor * sin_theta * np.sin(phi),
                0: prefactor * cos_theta,
                1: prefactor * sin_theta * np.cos(phi),
            }[m]
        if l == 2:
            return {
                -2: math.sqrt(15 / (16 * math.pi)) * sin_theta**2 * np.sin(2 * phi),
                -1: math.sqrt(15 / (4 * math.pi)) * cos_theta * sin_theta * np.sin(phi),
                0: 0.25 * math.sqrt(5 / math.pi) * (3 * cos_theta**2 - 1),
                1: math.sqrt(15 / (4 * math.pi)) * cos_theta * sin_theta * np.cos(phi),
                2: math.sqrt(15 / (16 * math.pi)) * sin_theta**2 * np.cos(2 * phi),
            }[m]
        if l == 3:
            return {
                -3: 0.25 * math.sqrt(35 / (2 * math.pi)) * sin_theta**3 * np.sin(3 * phi),
                -2: 0.25 * math.sqrt(105 / math.pi) * sin_theta**2 * cos_theta * np.sin(2 * phi),
                -1: 0.25 * math.sqrt(21 / (2 * math.pi)) * sin_theta * (5 * cos_theta**2 - 1) * np.sin(phi),
                0: 0.25 * math.sqrt(7 / math.pi) * (5 * cos_theta**3 - 3 * cos_theta),
                1: 0.25 * math.sqrt(21 / (2 * math.pi)) * sin_theta * (5 * cos_theta**2 - 1) * np.cos(phi),
                2: 0.25 * math.sqrt(105 / math.pi) * sin_theta**2 * cos_theta * np.cos(2 * phi),
                3: 0.25 * math.sqrt(35 / (2 * math.pi)) * sin_theta**3 * np.cos(3 * phi),
            }[m]
        raise ValueError(f"No real spherical harmonic for l={l}, m={m}.")
