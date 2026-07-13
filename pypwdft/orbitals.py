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

import numpy as np

def rotate_wavefunction_real(wavefunction):
    """
    Apply the global phase that makes a wavefunction maximally real.

    The returned complex array minimizes ``norm(result.imag)`` over all global
    phase rotations.  A global phase can make a Gamma-point wavefunction fully
    real when its imaginary part is only due to the arbitrary phase of an
    eigenvector.  It cannot remove a spatially varying phase or, in general, a
    complex mixture of **degenerate** orbitals.

    Args:
        wavefunction (array_like): Complex real- or reciprocal-space orbital.

    Returns:
        numpy.ndarray: A phase-rotated copy of ``wavefunction``.
    """
    psi = np.asarray(wavefunction, dtype=np.complex128)
    if psi.size == 0:
        return psi.copy()

    # For phi = exp(-i theta) psi, minimizing ||Im(phi)|| is equivalent to
    # maximizing Re(exp(-2i theta) sum(psi**2)).
    pseudo_norm = np.sum(psi * psi)
    if pseudo_norm == 0:
        return psi.copy()

    rotated = psi * np.exp(-0.5j * np.angle(pseudo_norm))

    # Resolve the remaining arbitrary sign reproducibly.
    pivot = np.unravel_index(np.argmax(np.abs(rotated)), rotated.shape)
    if rotated[pivot].real < 0:
        rotated = -rotated
    return rotated


def imaginary_fraction(wavefunction):
    """
    Return the fraction of a wavefunction's norm in its imaginary part.
    """
    psi = np.asarray(wavefunction)
    norm = np.linalg.norm(psi)
    if norm == 0:
        return 0.0
    return float(np.linalg.norm(psi.imag) / norm)
