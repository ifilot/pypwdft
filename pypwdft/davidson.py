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
"""

import numpy as np


def initial_guess(kinetic, nstates, xp=np):
    """
    Build a deterministic starting block: the ``nstates`` plane waves of
    lowest kinetic energy, plus a small fixed perturbation that prevents the
    block from being orthogonal to any eigenvector by symmetry.

    Args:
        kinetic: kinetic energy of every active plane wave, shape ``(npw,)``.
        nstates (int): number of states to solve for.
        xp (module, optional): array module, NumPy or CuPy.
    """
    kinetic = np.asarray(kinetic.get() if hasattr(kinetic, 'get') else kinetic)
    npw = len(kinetic)
    guess = np.random.default_rng(0).standard_normal((npw, nstates))
    guess *= 1e-2 / (1.0 + kinetic[:, None])
    lowest = np.argsort(kinetic, kind='stable')[:nstates]
    guess[lowest, np.arange(nstates)] += 1.0
    return xp.asarray(guess, dtype=np.complex128)


def davidson(apply, guess, kinetic, tol, maxiter=200, xp=np):
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

    Returns:
        tuple: eigenvalues ``(nstates,)`` in ascending order, orthonormal
        eigenvectors ``(npw, nstates)`` and whether all residuals converged.
    """
    nstates = guess.shape[1]
    max_space = max(4 * nstates, nstates + 16)
    kinetic = kinetic[:, None]

    V, _ = xp.linalg.qr(guess)      # orthonormal basis of the search space
    HV = apply(V)
    converged = False

    for _ in range(maxiter):
        # Rayleigh-Ritz: the best eigenpair estimates within the search space
        ritz_values, ritz_vectors = xp.linalg.eigh(V.conj().T @ HV)
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
        T = residual[:, active] * (poly / (poly + 16 * x**4))

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
            T -= V @ (V.conj().T @ T)
        T, r = xp.linalg.qr(T)
        T = T[:, xp.abs(xp.diagonal(r)) > 1e-8]
        if T.shape[1] == 0:
            break

        V = xp.concatenate([V, T], axis=1)
        HV = xp.concatenate([HV, apply(T)], axis=1)

    return theta, X, converged
