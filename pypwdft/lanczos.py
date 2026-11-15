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
Lowest eigenpairs with a Lanczos solver (ARPACK or CuPy ``eigsh``), made
reliable for degenerate spectra.

A Lanczos solver builds its search space from a single start vector, which has
two consequences that a plain call to ``eigsh`` does not guard against:

* A start vector with a spatial symmetry never finds states of another
  symmetry. A single plane wave, for instance, has no component along the pi
  orbitals of a linear molecule.
* The search space contains only *one* vector of every degenerate level. The
  other members appear only through rounding errors, so at any finite
  tolerance the solver may return a higher state in their place, without any
  warning.

:func:`lanczos_lowest` therefore starts from a generic vector and verifies
its result by deflation: the states found are shifted out of the way and the
solver is asked for the lowest remaining state. Anything that turns up below
the highest state found was missed, and is added. The iteration count of every
solve is bounded, so that a stagnating solve (which happens when the requested
number of states cuts through a degenerate level) is detected rather than
endured; it is then repeated with one more state.
"""

import numpy as np

from .fft import limit_blas_threads

# Maximum number of restarts of one Lanczos solve. Converging solves need
# far fewer; the bound exists to detect stagnation.
MAX_RESTARTS = 100

# Maximum number of extra states requested to get past a stagnating solve.
MAX_EXTRA_STATES = 4

# Loosest tolerances used, whatever the caller asks for. A Lanczos solve
# stops as soon as the states it has are converged; with a loose tolerance
# that can happen before a lower state has emerged in the search space. The
# verification solves ask for a single state and are cheap, so they are run
# with a tighter tolerance.
MAX_TOL = 1e-4
MAX_VERIFY_TOL = 1e-6


def start_vector(kinetic, xp=np):
    """
    Deterministic start vector without any spatial symmetry.

    Random complex coefficients, damped at high kinetic energy so that the
    vector is dominated by the smooth functions the lowest states consist of.

    Args:
        kinetic: kinetic energy of every active plane wave, shape ``(npw,)``.
        xp (module, optional): array module, NumPy or CuPy.
    """
    kinetic = np.asarray(kinetic.get() if hasattr(kinetic, 'get') else kinetic)
    rng = np.random.default_rng(0)
    vector = rng.standard_normal(len(kinetic)) \
           + 1j * rng.standard_normal(len(kinetic))
    vector /= 1.0 + kinetic
    return xp.asarray(vector / np.linalg.norm(vector))


def lanczos_lowest(sparse_linalg, apply, nstates, v0, tol=0, xp=np,
                   partner=None):
    """
    Find the ``nstates`` lowest eigenpairs of a Hermitian operator.

    Args:
        sparse_linalg (module): ``scipy.sparse.linalg`` or
            ``cupyx.scipy.sparse.linalg``; provides ``eigsh``.
        apply (callable): applies the operator to a block of column vectors
            of shape ``(n, m)``.
        nstates (int): number of eigenpairs.
        v0: start vector of shape ``(n,)``. It must not have a spatial
            symmetry, see :func:`start_vector`.
        tol (float, optional): tolerance passed to ``eigsh``; zero requests
            machine precision.
        xp (module, optional): array module owning the arrays.
        partner (optional): index of ``-G`` for every plane wave, see
            :func:`pypwdft.davidson.conjugation_partner`. When given, the
            operator must be real in real space; the eigenvectors returned
            are then real functions in real space, also within degenerate
            levels and for a level that ``nstates`` cuts through.

    Returns:
        tuple: eigenvalues ``(nstates,)`` in ascending order and orthonormal
        eigenvectors ``(n, nstates)``.

    Raises:
        RuntimeError: if the solver does not converge.
    """
    if xp is not np:
        return _lanczos_lowest(
            sparse_linalg, apply, nstates, v0, tol, xp, partner
        )
    # ARPACK calls back into NumPy for every operator application. SciPy and
    # NumPy each bring their own BLAS thread pool, and two pools taking turns
    # on small products slow each other down by an order of magnitude, so
    # BLAS is kept to a single thread for the duration of the solve.
    with limit_blas_threads(2):
        return _lanczos_lowest(
            sparse_linalg, apply, nstates, v0, tol, xp, partner
        )


def _lanczos_lowest(sparse_linalg, apply, nstates, v0, tol, xp, partner):
    """Implementation of :func:`lanczos_lowest`."""
    size = len(v0)
    tol = min(tol, MAX_TOL)
    margin = max(tol, 1e-9)

    def solve(operator, count, start, tol=tol):
        """One bounded ``eigsh`` call; ``None`` if it did not converge."""
        def matvec(vector):
            return operator(vector.reshape(-1, 1))[:, 0]

        linear_operator = sparse_linalg.LinearOperator(
            (size, size), matvec=matvec, dtype=np.complex128
        )
        # eigsh may overwrite its start vector, so it is given a copy
        start = xp.ascontiguousarray(start).copy()
        if xp is np:
            try:
                values, vectors = sparse_linalg.eigsh(
                    linear_operator, count, which='SA', v0=start, tol=tol,
                    maxiter=MAX_RESTARTS,
                )
            except sparse_linalg.ArpackNoConvergence:
                return None
        else:
            # CuPy counts Lanczos steps rather than restarts and returns its
            # last iterate without complaint, so convergence is checked here.
            ncv = min(max(2 * count, count + 32), size - 1)
            values, vectors = sparse_linalg.eigsh(
                linear_operator, count, which='SA', v0=start, tol=tol,
                maxiter=MAX_RESTARTS * (ncv - count),
            )
            residual = xp.linalg.norm(
                operator(vectors) - vectors * values, axis=0
            )
            scale = max(1.0, float(xp.abs(values).max()))
            if float(residual.max()) > 100 * max(tol, 1e-10) * scale:
                return None
        return values, vectors

    # Solve, asking for one more state whenever the solve stagnates.
    for extra in range(MAX_EXTRA_STATES + 1):
        count = nstates + extra
        if count >= size - 1:
            break
        result = solve(apply, count, v0)
        if result is not None:
            break
    else:
        result = None
    if result is None:
        raise RuntimeError(
            'The Lanczos eigensolver did not converge. '
            "Use eigensolver='davidson' instead."
        )
    values, vectors = result

    # Verify by deflation: shift the states found above the highest of them
    # and look for the lowest state that remains. Repeat while that state
    # lies below the highest wanted state, i.e. while a state was missed.
    probe = start_vector_like(v0, xp)
    for _ in range(nstates + MAX_EXTRA_STATES + 1):
        order = xp.argsort(values)
        highest = float(values[order[nstates - 1]])
        shift = float(values.max() - values.min()) + 1.0
        found = vectors

        def deflated(block, found=found, shift=shift):
            return apply(block) + shift * (found @ (found.conj().T @ block))

        result = solve(
            deflated, 1, probe - found @ (found.conj().T @ probe),
            tol=min(tol, MAX_VERIFY_TOL),
        )
        if result is None:
            raise RuntimeError(
                'The Lanczos eigensolver did not converge while verifying '
                "its result. Use eigensolver='davidson' instead."
            )
        value, vector = result
        if float(value[0]) >= highest - margin * max(1.0, abs(highest)):
            break
        # A missed state: orthogonalise it against the others and keep it.
        vector = vector - found @ (found.conj().T @ vector)
        vector /= xp.linalg.norm(vector)
        values = xp.concatenate([values, value])
        vectors = xp.concatenate([vectors, vector], axis=1)
    else:
        raise RuntimeError(
            'The Lanczos eigensolver keeps finding missed states. '
            "Use eigensolver='davidson' instead."
        )

    # Rayleigh-Ritz within everything found, then keep the lowest states.
    # This is done unconditionally: for a complex Hermitian operator ARPACK
    # runs its general (Arnoldi) algorithm, whose eigenvectors within a
    # degenerate level are not orthogonal to each other.
    if partner is None:
        basis, _ = xp.linalg.qr(vectors)
        projected = basis.conj().T @ apply(basis)
        projected = 0.5 * (projected + projected.conj().T)
    else:
        basis = _real_basis(vectors, partner, xp)
        # between real functions the operator is real symmetric
        projected = (basis.conj().T @ apply(basis)).real
        projected = 0.5 * (projected + projected.T)
    values, rotation = xp.linalg.eigh(projected)
    return values[:nstates], basis @ rotation[:, :nstates]


def _real_basis(vectors, partner, xp):
    """
    Orthonormal basis of real-space-real functions spanning ``vectors``.

    With ``K`` the complex conjugation in real space, ``(K c)(G) =
    conj(c(-G))``, an eigenvector ``v`` of a real operator has the real
    partners ``(v + K v) / 2`` and ``i (v - K v) / 2``, which are eigenvectors
    of the same level and together reproduce ``v``. For a non-degenerate
    state the two are parallel; for a complex mixture within a degenerate
    level they are independent. The dependent directions are removed.
    """
    conjugate = vectors[partner].conj()
    candidates = xp.concatenate(
        [0.5 * (vectors + conjugate), 0.5j * (vectors - conjugate)], axis=1
    )
    # The overlap of real functions is real, so the orthonormalisation only
    # takes real combinations and the basis stays real in real space.
    overlap = (candidates.conj().T @ candidates).real
    weights, rotation = xp.linalg.eigh(0.5 * (overlap + overlap.T))
    # Keep the significant directions, but never fewer than the number of
    # vectors: their span has at least that dimension.
    keep = weights > 1e-4 * float(weights[-1])
    keep[-vectors.shape[1]:] = True
    return (candidates @ rotation[:, keep]) / xp.sqrt(weights[keep])


def start_vector_like(vector, xp=np):
    """
    Deterministic generic vector with the shape of ``vector``, used as the
    start of the verification solves.
    """
    rng = np.random.default_rng(1)
    probe = rng.standard_normal(len(vector)) \
          + 1j * rng.standard_normal(len(vector))
    return xp.asarray(probe / np.linalg.norm(probe))
