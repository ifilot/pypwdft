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
Density-mixing algorithms used by the self-consistent field solver.

The Kohn--Sham equations define a fixed-point problem for the electron
density.  Starting from an input density, the Hamiltonian produces orbitals
and therefore a new output density.  Directly using that output as the next
input often causes oscillation.  This module provides damped linear mixing and
Pulay mixing, also known as direct inversion in the iterative subspace (DIIS).
"""

from __future__ import annotations

import numpy as np


class DensityMixer:
    r"""
    Mix SCF densities using linear or Pulay (DIIS) extrapolation.

    For input density :math:`\rho_i^\mathrm{in}` and output density
    :math:`\rho_i^\mathrm{out}`, the fixed-point residual is

    .. math::

        R_i = \rho_i^\mathrm{out} - \rho_i^\mathrm{in}.

    Linear mixing forms the next input as

    .. math::

        \rho_{i+1}^\mathrm{in} = \rho_i^\mathrm{in} + \beta R_i,

    where :math:`\beta` is ``fraction``. Pulay mixing retains several recent
    residuals and finds the normalized linear combination having the smallest
    residual norm. The first iteration and unsafe or ill-conditioned Pulay
    extrapolations automatically fall back to linear mixing.

    Args:
        xp: Array module used by the solver. This is normally :mod:`numpy`, or
            :mod:`cupy` for GPU calculations.
        method (str, optional): ``"pulay"`` or ``"linear"``. Defaults to
            ``"pulay"``.
        fraction (float, optional): Linear damping parameter :math:`\beta` in
            the interval ``(0, 1]``. Defaults to 0.5.
        history (int, optional): Maximum number of density-residual pairs kept
            in the Pulay subspace. Defaults to 6.
        target_mean (float, optional): Required spatial mean of the density.
            In the SCF solver this equals the electron count divided by the
            cell volume. When supplied, every returned density is normalized
            to this value.
    """

    def __init__(
        self,
        xp,
        *,
        method="pulay",
        fraction=0.5,
        history=6,
        target_mean=None,
    ):
        # Normalize and validate all user-facing settings before allocating
        # history.  A Pulay subspace needs at least two residual vectors.
        method = str(method).lower()
        if method not in {"linear", "pulay"}:
            raise ValueError("mixing must be 'linear' or 'pulay'.")
        if not 0 < fraction <= 1:
            raise ValueError("mixing_fraction must be in the interval (0, 1].")
        if not isinstance(history, (int, np.integer)) or history < 2:
            raise ValueError("mixing_history must be an integer of at least 2.")

        # Keep the numerical array module rather than importing CuPy here. This
        # lets the same implementation operate on CPU and GPU density arrays.
        self.xp = xp
        self.method = method
        self.fraction = float(fraction)
        self.history = int(history)
        self.target_mean = target_mean

        # The two lists form a synchronized, first-in-first-out history. Each
        # density at index i belongs to the residual at the same index.
        self._densities = []
        self._residuals = []

        # Expose which algorithm supplied the most recent update. This is also
        # useful for testing whether a Pulay proposal passed the safeguards.
        self.last_step = "linear"

    def update(self, input_density, output_density):
        """
        Construct the density used as input to the next SCF iteration.

        Args:
            input_density (array_like): Density used to construct the current
                Kohn--Sham Hamiltonian.
            output_density (array_like): Density constructed from the current
                occupied Kohn--Sham orbitals.

        Returns:
            array_like: Mixed and charge-normalized next input density. The
            returned array belongs to the same numerical backend as the input.
        """
        xp = self.xp

        # The residual vanishes at self-consistency. Linear damping is always
        # calculated because it is both the first step and the safe fallback.
        residual = output_density - input_density
        damped_density = input_density + self.fraction * residual

        # Store copies because the SCF driver reuses its density arrays. Keep a
        # bounded history so memory and the small DIIS solve remain predictable.
        self._densities.append(damped_density.copy())
        self._residuals.append(residual.copy())
        if len(self._densities) > self.history:
            self._densities.pop(0)
            self._residuals.pop(0)

        # Linear mixing is either explicitly requested or required until two
        # independent residuals are available to define a Pulay subspace.
        if self.method == "linear" or len(self._residuals) < 2:
            self.last_step = "linear"
            return self._normalise(damped_density)

        # Find the combination of history vectors that minimizes the residual
        # norm while preserving a coefficient sum of one. A failed solve is
        # deliberately non-fatal: the damped linear step remains well defined.
        coefficients = self._pulay_coefficients()
        if coefficients is None:
            self.last_step = "linear"
            return self._normalise(damped_density)

        # Apply the DIIS coefficients to the damped densities. The coefficients
        # are ordinary CPU scalars, whereas the density may reside on a GPU.
        candidate = xp.zeros_like(input_density)
        for coefficient, density in zip(coefficients, self._densities):
            candidate += coefficient * density

        # A nearly singular history can yield a mathematically valid solution
        # with very large cancelling coefficients. Reject a Pulay proposal if
        # it contains NaN/Inf values, becomes negative, or moves more than ten
        # times farther than the corresponding damped linear step. Negative
        # densities are especially problematic for fractional powers in LDA.
        candidate_step = self._rms(candidate - input_density)
        linear_step = self._rms(damped_density - input_density)
        safe_step = candidate_step <= max(10.0 * linear_step, 1e-14)
        finite = self._boolean(xp.all(xp.isfinite(candidate)))
        nonnegative = self._scalar(xp.min(candidate)) >= 0.0
        if not finite or not nonnegative or not safe_step:
            self.last_step = "linear"
            return self._normalise(damped_density)

        self.last_step = "pulay"
        return self._normalise(candidate)

    def _pulay_coefficients(self):
        r"""
        Solve the constrained residual-minimization problem on the CPU.

        If :math:`B_{ij}=\langle R_i,R_j\rangle`, minimizing the residual norm
        subject to :math:`\sum_i c_i=1` gives the augmented linear system

        .. math::

            \begin{pmatrix} B & 1 \\ 1^T & 0 \end{pmatrix}
            \begin{pmatrix} c \\ \lambda \end{pmatrix}
            = \begin{pmatrix} 0 \\ 1 \end{pmatrix}.

        Only this small history matrix is solved on the CPU. For a GPU
        calculation, forming its entries transfers only scalar inner products,
        not full density grids.

        Returns:
            numpy.ndarray or None: Pulay coefficients, or ``None`` when the
            history cannot provide a safe extrapolation.
        """
        count = len(self._residuals)

        # Construct the symmetric residual-overlap matrix. The spatial mean is
        # used rather than a sum; this changes only a common scale and keeps
        # matrix magnitudes comparable for different FFT-grid sizes.
        overlap = np.empty((count, count), dtype=float)
        for i, residual_i in enumerate(self._residuals):
            for j in range(i + 1):
                value = self._scalar(
                    self.xp.mean(residual_i * self._residuals[j])
                )
                overlap[i, j] = overlap[j, i] = value

        # Scale by the largest residual norm before adding a small diagonal
        # regularizer. The regularizer makes nearly dependent histories less
        # singular without materially changing a well-conditioned subspace.
        scale = float(np.max(np.abs(np.diag(overlap))))
        if not np.isfinite(scale) or scale <= np.finfo(float).tiny:
            return None
        overlap /= scale
        overlap.flat[::count + 1] += 1e-10

        # Add a Lagrange multiplier as the final row and column to enforce that
        # the density-combination coefficients sum exactly to one.
        system = np.empty((count + 1, count + 1), dtype=float)
        system[:count, :count] = overlap
        system[:count, count] = 1.0
        system[count, :count] = 1.0
        system[count, count] = 0.0
        rhs = np.zeros(count + 1, dtype=float)
        rhs[count] = 1.0
        # DIIS acceleration is optional. If numerical linear algebra detects a
        # singular system, signal the caller to use the linear fallback.
        try:
            coefficients = np.linalg.solve(system, rhs)[:count]
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(coefficients)):
            return None
        # Very large individual coefficients indicate cancellation between
        # almost redundant history vectors and tend to produce unstable steps.
        if np.max(np.abs(coefficients)) > 20.0:
            return None
        return coefficients

    def _normalise(self, density):
        """
        Restore the prescribed electron count after numerical extrapolation.

        Pulay coefficients sum to one, so charge should be conserved already.
        Explicit normalization removes small floating-point drift, which is
        particularly useful over long SCF cycles.
        """
        if self.target_mean is None:
            return density
        mean = self._scalar(self.xp.mean(density))
        if not np.isfinite(mean) or mean <= 0:
            raise RuntimeError("density mixing produced an invalid electron density.")
        return density * (self.target_mean / mean)

    def _rms(self, value):
        """Return the root-mean-square magnitude as a Python float."""
        return self._scalar(self.xp.sqrt(self.xp.mean(value * value)))

    @staticmethod
    def _scalar(value):
        """Convert a NumPy or CuPy scalar to an ordinary Python float."""
        # CuPy scalars provide ``get`` to copy their value from device to host;
        # NumPy scalars do not and can be converted directly.
        if hasattr(value, "get"):
            value = value.get()
        return float(np.real(value))

    @staticmethod
    def _boolean(value):
        """Convert a NumPy or CuPy Boolean scalar to a Python bool."""
        if hasattr(value, "get"):
            value = value.get()
        return bool(value)
