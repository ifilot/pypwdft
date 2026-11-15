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
    complex mixture of **degenerate** orbitals; use
    :func:`rotate_orbitals_real` for a set of orbitals with their energies.

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


def degenerate_groups(energies, tol=1e-4):
    """
    Group orbital indices whose energies lie within ``tol`` of a neighbour.

    Args:
        energies (array_like): Orbital energies in Hartree.
        tol (float, optional): Largest energy gap inside one group.

    Returns:
        list of numpy.ndarray: Index arrays, ordered by increasing energy.
    """
    energies = np.asarray(energies, dtype=float)
    if energies.size == 0:
        return []
    order = np.argsort(energies, kind='stable')
    breaks = np.flatnonzero(np.diff(energies[order]) > tol) + 1
    return np.split(order, breaks)


def rotate_orbitals_real(orbitals, energies=None, degeneracy_tol=1e-4):
    """
    Make a set of orbitals maximally real, including degenerate ones.

    A non-degenerate Gamma-point orbital is real up to a global phase.  Within
    a degenerate level, however, an eigensolver may return any unitary mixture
    of the real eigenfunctions, for example ``(p_x + i p_y) / sqrt(2)``, and no
    phase factor makes such an orbital real.  This function therefore treats
    every degenerate level as a whole: it applies the unitary transformation
    within the level that minimizes the summed squared norm of the imaginary
    parts.  The transformed orbitals span the same level, so they remain
    eigenfunctions with the same energy, and they remain orthonormal.

    For a level of ``m`` orbitals ``psi_i`` the transformation follows from
    the Takagi factorization of the complex-symmetric matrix
    ``M_ij = sum(psi_i * psi_j)`` (no complex conjugation).  For ``m = 1``
    this reduces to :func:`rotate_wavefunction_real`.

    A level can only be made fully real if all of its members are present.
    If the set of orbitals cuts through a degenerate level, the members that
    are present generally keep an imaginary part; :func:`imaginary_fraction`
    reports it.

    Args:
        orbitals (array_like): Orbitals with shape ``(norb, ...)``, in real
            space.
        energies (array_like, optional): Orbital energies in Hartree, used to
            identify degenerate levels.  Without energies every orbital is
            treated as non-degenerate.
        degeneracy_tol (float, optional): Largest energy gap, in Hartree,
            between neighbouring members of one level.  Orbitals within this
            tolerance that are already real up to a phase are left unmixed.

    Returns:
        numpy.ndarray: Transformed copy of ``orbitals``.
    """
    psi = np.array(orbitals, dtype=np.complex128)
    if psi.ndim < 1:
        raise ValueError("orbitals must have shape (norb, ...).")
    if energies is None:
        groups = [np.array([index]) for index in range(len(psi))]
    else:
        if len(energies) != len(psi):
            raise ValueError("energies must contain one value per orbital.")
        groups = degenerate_groups(energies, degeneracy_tol)

    flat = psi.reshape(len(psi), -1)
    for group in groups:
        block = flat[group]
        if block.size == 0:
            continue
        block = _takagi_rotation(block @ block.T).T @ block
        # Resolve the remaining arbitrary sign of every orbital reproducibly.
        pivots = np.argmax(np.abs(block), axis=1)
        signs = np.sign(block[np.arange(len(block)), pivots].real)
        signs[signs == 0] = 1.0
        flat[group] = block * signs[:, None]
    return psi


def _takagi_rotation(pseudo_overlap):
    """
    Return the unitary ``W`` for which ``psi @ W`` is maximally real.

    ``pseudo_overlap`` is the complex-symmetric matrix ``M = psi^T psi`` of
    the orbitals stored as columns.  Minimizing the imaginary norm of the
    rotated orbitals is equivalent to maximizing ``Re tr(W^T M W)``, which is
    solved by the Takagi factorization ``M = V S V^T``: ``W = conj(V)``.
    """
    matrix = 0.5 * (pseudo_overlap + pseudo_overlap.T)
    size = len(matrix)
    scale = np.linalg.norm(matrix)
    if scale == 0:
        return np.eye(size, dtype=np.complex128)

    # Takagi vectors v = x + iy satisfy M conj(v) = s v, which is the real
    # symmetric eigenproblem [[A, B], [B, -A]] (x, y) = s (x, y) for
    # M = A + iB.  Its eigenvalues come in pairs +s and -s; the Takagi
    # vectors belong to the non-negative ones.  Unlike a route via the SVD,
    # this stays well defined when singular values coincide, which is the
    # normal situation here: a complete degenerate level has M unitary.
    real, imag = matrix.real, matrix.imag
    _, vectors = np.linalg.eigh(np.block([[real, imag], [imag, -real]]))
    top = vectors[:, size:]
    takagi = top[:size] + 1j * top[size:]
    # Vectors of (near-)zero singular values need not be orthogonal to the
    # others; the polar factor restores exact unitarity and is a no-op
    # otherwise.
    left, _, right = np.linalg.svd(takagi)
    rotation = (left @ right).conj()

    # Any real orthogonal matrix may still be applied.  Choose the one that
    # keeps the orbitals closest to the originals, each with its own best
    # phase, so that orbitals which are already real up to a phase are not
    # mixed with their (near-)degenerate partners.
    phases = np.exp(-0.5j * np.angle(np.diagonal(matrix)))
    left, _, right = np.linalg.svd((rotation.conj().T * phases).real)
    return rotation @ (left @ right)


def imaginary_fraction(wavefunction):
    """
    Return the fraction of a wavefunction's norm in its imaginary part.
    """
    psi = np.asarray(wavefunction)
    norm = np.linalg.norm(psi)
    if norm == 0:
        return 0.0
    return float(np.linalg.norm(psi.imag) / norm)
