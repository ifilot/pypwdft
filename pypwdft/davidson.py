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

"""
Block Davidson eigensolver for the lowest states of a Hermitian operator.

Unlike a Lanczos solver, Davidson starts from a full block of trial vectors,
so the orbitals of the previous SCF iteration can be reused and only a few
Hamiltonian applications per band are needed.

Real orbitals
-------------
At the Gamma point the Kohn-Sham Hamiltonian is real in real space, so its
eigenfunctions can be chosen real. In plane-wave coefficients ``c(G)`` a real
function obeys ``c(-G) = conj(c(G))``. Given the index of ``-G`` for every
plane wave (see :func:`conjugation_partner`), the solver works entirely within
this set of vectors: the Hamiltonian maps it onto itself, the projected
eigenproblem becomes real symmetric, and every returned eigenvector is a real
function in real space. This also holds for each individual member of a
degenerate level, which an unconstrained solver may return as an arbitrary
complex mixture.
"""

import numpy as np


def conjugation_partner(pw_mask):
    """
    For every active plane wave ``G``, find the active index of ``-G``.

    Args:
        pw_mask (array_like): boolean mask of the active plane waves on the
            cubic FFT grid, in FFT (``fftfreq``) ordering.

    Returns:
        numpy.ndarray or None: integer array ``partner`` of shape ``(npw,)``
        such that plane wave ``partner[i]`` is the negative of plane wave
        ``i``; ``None`` if the active set is not closed under ``G -> -G``.
    """
    mask = np.asarray(pw_mask, dtype=bool)
    active = np.flatnonzero(mask)
    position = np.full(mask.size, -1, dtype=np.int64)
    position[active] = np.arange(len(active))
    # On an FFT grid the index of -G along one axis is (-i) mod n.
    indices = np.unravel_index(active, mask.shape)
    negated = np.ravel_multi_index(
        tuple((-index) % size for index, size in zip(indices, mask.shape)),
        mask.shape,
    )
    partner = position[negated]
    if np.any(partner < 0):
        return None
    return partner


def initial_guess(kinetic, nstates, xp=np, partner=None):
    """
    Build a deterministic starting block: the ``nstates`` plane waves of
    lowest kinetic energy, plus a small fixed perturbation that prevents the
    block from being orthogonal to any eigenvector by symmetry.

    Args:
        kinetic: kinetic energy of every active plane wave, shape ``(npw,)``.
        nstates (int): number of states to solve for.
        xp (module, optional): array module, NumPy or CuPy.
        partner (array_like, optional): index of ``-G`` for every plane wave,
            see :func:`conjugation_partner`. When given, the block consists
            of real functions: the cosine and sine combinations of the
            lowest plane-wave pairs ``(G, -G)``.
    """
    kinetic = np.asarray(kinetic.get() if hasattr(kinetic, 'get') else kinetic)
    npw = len(kinetic)
    guess = np.random.default_rng(0).standard_normal((npw, nstates))
    guess *= 1e-2 / (1.0 + kinetic[:, None])
    order = np.argsort(kinetic, kind='stable')
    if partner is None:
        guess[order[:nstates], np.arange(nstates)] += 1.0
        return xp.asarray(guess, dtype=np.complex128)

    partner = np.asarray(partner.get() if hasattr(partner, 'get') else partner)
    guess = 0.5 * (guess + guess[partner]).astype(np.complex128)
    # The same plane waves as above, in the same order, each replaced by a
    # real function: cos(G.r) the first time a pair (G, -G) is met and
    # sin(G.r) the second time. Complete pairs thus span the same space as
    # the two plane waves themselves.
    met = np.zeros(npw, dtype=bool)
    for column, index in enumerate(order[:nstates]):
        other = partner[index]
        if other == index:                      # G = 0
            guess[index, column] += 1.0
        elif not met[other]:                    # cos(G.r)
            guess[index, column] += np.sqrt(0.5)
            guess[other, column] += np.sqrt(0.5)
        else:                                   # sin(G.r)
            guess[index, column] += 1j * np.sqrt(0.5)
            guess[other, column] -= 1j * np.sqrt(0.5)
        met[index] = True
    return xp.asarray(guess)


def davidson(apply, guess, kinetic, tol, maxiter=200, xp=np, partner=None):
    """
    Find the lowest eigenpairs of a Hermitian operator.

    Args:
        apply (callable): applies the operator to every column of a block
            of shape ``(npw, n)``.
        guess: starting block of shape ``(npw, nstates)``; its column count
            sets the number of eigenpairs.
        kinetic: kinetic energy of every active plane wave, shape ``(npw,)``,
            used by the preconditioner.
        tol (float): convergence threshold on the residual norm of each
            eigenpair.
        maxiter (int, optional): maximum number of subspace expansions.
        xp (module, optional): array module owning the arrays.
        partner (optional): index of ``-G`` for every plane wave, see
            :func:`conjugation_partner`. When given, the operator must be
            real in real space; the eigenvectors are then real functions in
            real space, also within degenerate levels.

    Returns:
        tuple: eigenvalues ``(nstates,)`` in ascending order, orthonormal
        eigenvectors ``(npw, nstates)`` and whether all residuals converged.
    """
    nstates = guess.shape[1]
    max_space = max(4 * nstates, nstates + 16)
    kinetic = kinetic[:, None]

    def real_part(block):
        """Project a block onto the vectors with c(-G) = conj(c(G))."""
        if partner is None:
            return block
        return 0.5 * (block + block[partner].conj())

    # orthonormal basis of the search space; QR keeps real functions real
    V, _ = xp.linalg.qr(real_part(guess))
    HV = apply(V)
    converged = False

    for _ in range(maxiter):
        # Rayleigh-Ritz: the best eigenpair estimates within the search space
        projected = V.conj().T @ HV
        if partner is not None:
            # Between real functions the Hamiltonian matrix is real symmetric;
            # dropping the rounding-level imaginary part makes the Ritz
            # vectors real combinations of real functions.
            projected = 0.5 * (projected.real + projected.real.T)
        ritz_values, ritz_vectors = xp.linalg.eigh(projected)
        theta, y = ritz_values[:nstates], ritz_vectors[:, :nstates]
        X = V @ y
        HX = HV @ y

        residual = HX - X * theta
        active = xp.linalg.norm(residual, axis=0) > tol
        if not bool(active.any()):
            converged = True
            break

        # Teter-Payne-Allan preconditioner: damps the high-kinetic-energy
        # components that dominate the residual of a plane-wave Hamiltonian.
        ekin = xp.real(xp.sum(xp.abs(X[:, active])**2 * kinetic, axis=0))
        x = kinetic / (1.5 * ekin)
        poly = 27 + 18 * x + 12 * x**2 + 8 * x**3
        T = real_part(residual[:, active] * (poly / (poly + 16 * x**4)))

        # Thick restart: when the search space is full, collapse it onto its
        # lowest Ritz vectors. Keeping twice the number of wanted states
        # preserves the directions that the slowest bands still need.
        if V.shape[1] + T.shape[1] > max_space:
            keep = ritz_vectors[:, :2 * nstates]
            V, HV = V @ keep, HV @ keep

        # Orthogonalise the corrections against the search space (twice for
        # numerical stability) and among themselves; drop dependent ones.
        T /= xp.linalg.norm(T, axis=0)
        for _ in range(2):
            overlap = V.conj().T @ T
            if partner is not None:
                overlap = overlap.real
            T -= V @ overlap
        T, r = xp.linalg.qr(T)
        T = T[:, xp.abs(xp.diagonal(r)) > 1e-8]
        if T.shape[1] == 0:
            break

        V = xp.concatenate([V, T], axis=1)
        HV = xp.concatenate([HV, apply(T)], axis=1)

    return theta, X, converged
