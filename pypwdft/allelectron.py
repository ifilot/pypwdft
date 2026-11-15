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

"""All-electron ionic model: bare nuclear Coulomb potentials."""

import numpy as np

from .elements import symbol


class AllElectron:
    """Bare-nucleus ionic model associated with a :class:`PeriodicSystem`.

    Exposes the same interface as :class:`GTHPseudopotential`, so the solver
    does not need to distinguish between the two ionic models.
    """

    label = "all-electron"

    def __init__(self, system):
        self.system = system

    @property
    def symbols(self):
        """Element symbols in atom order."""
        return tuple(symbol(z) for z in self.system.get_atom_charges())

    @property
    def nelec(self):
        """Number of electrons, equal to the total nuclear charge."""
        return self.system.get_nelec()

    def local_potential(self):
        """Return the nuclear attraction potential on the real-space FFT grid."""
        kvec = self.system.get_pw_k()
        k2 = self.system.get_pw_k2()

        # calculate structure factor
        sf = np.exp(1j * kvec @ self.system.get_atom_positions().T)

        with np.errstate(divide='ignore', invalid='ignore'):
            nucpotg = -4.0 * np.pi / k2
            nucpotg[0,0,0] = 0

        # produce the nuclear attraction field
        return np.fft.fftn(np.einsum('ijk,ijkl,l->ijk',
                                     nucpotg,
                                     sf,
                                     self.system.get_atom_charges())
                           ) / self.system.get_omega()

    def nonlocal_operator(self, array_module=np):
        """Return ``None``: the bare Coulomb potential is purely local."""
        return None

    def ionic_energy(self):
        """Return the periodic nucleus-nucleus repulsion energy."""
        return self.system.calculate_ewald_sum()
