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
Exchange-correlation functionals for spin-unpolarized densities.

Every functional returns the XC energy per particle and the multiplicative XC
potential on the real-space grid. ``xp`` is the array module (NumPy or CuPy)
that owns the density.
"""

import numpy as np

# Supported functionals and the aliases accepted for them.
FUNCTIONALS = ('lda', 'pbe')
_ALIASES = {'svwn5': 'lda'}


def normalise(functional:str) -> str:
    """
    Return the canonical name of a functional, resolving aliases.
    """
    name = str(functional).lower()
    name = _ALIASES.get(name, name)
    if name not in FUNCTIONALS:
        raise ValueError(
            f'Unknown exchange-correlation functional: {functional}'
        )
    return name


def evaluate(functional:str, rho:np.ndarray, kvec:np.ndarray=None, xp=np):
    """
    Evaluate XC energy per particle and its multiplicative potential.

    Args:
        functional (str): canonical functional name, see :func:`normalise`.
        rho (np.ndarray): electron density on the real-space grid.
        kvec (np.ndarray, optional): plane-wave vectors of the FFT grid,
            required by gradient-corrected functionals.
        xp (module, optional): array module owning ``rho`` and ``kvec``.
    """
    if functional == 'lda':
        ex, vx = lda_x(rho)
        ec, vc = lda_c_vwn(rho, xp)
        return ex + ec, vx + vc
    return pbe(rho, kvec, xp)


def lda_x(rho:np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Slater exchange functional, see Parr and Yang pages 154 and 155,
    equations 7.4.5 and 7.4.9
    """
    f = -3 / 4 * (3 / (2 * np.pi))**(2 / 3)
    rs = (3 / (4 * np.pi * rho))**(1 / 3)

    ex = f / rs
    vx = 4 / 3 * ex

    return ex, vx


def lda_c_vwn(rho:np.ndarray, xp=np) -> (np.ndarray, np.ndarray):
    """
    Vosko-Wilk-Nusair correlation functional, see Parr and Yang page 275
    equation E.27
    """
    # VWN5 fit to the Ceperley-Alder Monte Carlo electron-gas data.
    # A is expressed in Rydberg here and the factor 1/2 below converts
    # the correlation energy to Hartree.
    A = 0.0621814
    x0 = -0.10498
    b = 3.72744
    c = 12.9352

    rs = (3 / (4 * np.pi * rho))**(1 / 3)

    x = rs**(1/2)
    X = x**2 + b * x + c
    X0 = x0**2 + b * x0 + c
    fx0 = b * x0 / (x0**2 + b * x0 + c)
    tx = 2 * x + b
    Q = (4 * c - b**2)**(1/2)
    atan = xp.arctan(Q / (2*x+b))

    ec = A/2 * (xp.log(x**2/X) + 2*b/Q * atan - b*x0/X0 * \
                (xp.log((x-x0)**2 / X) + 2 * (b + 2*x0) / Q * atan))

    tt = tx**2 + Q**2
    vc = ec - x * A / 12 * (2 / x - tx / X - 4 * b / tt - fx0 * \
        (2 / (x - x0) - tx / X - 4 * (2 * x0 + b) / tt))

    return ec,vc


def pbe(rho:np.ndarray, kvec:np.ndarray, xp=np) -> (np.ndarray, np.ndarray):
    """
    Spin-unpolarized PBE energy per particle and GGA potential.
    """
    rho = xp.maximum(xp.asarray(rho, dtype=float), 1e-14)
    grad = _gradient(rho, kvec, xp)
    norm_grad = xp.linalg.norm(grad, axis=-1)

    # PBE exchange: LDA exchange times the enhancement factor F_x(s).
    kappa = 0.804
    mu = 0.2195149727645171
    kf = (3 * np.pi**2 * rho)**(1 / 3)
    s = norm_grad / (2 * kf * rho)
    denominator = 1 + mu * s**2 / kappa
    enhancement_correction = kappa - kappa / denominator
    ex_lda = -3 * kf / (4 * np.pi)
    ex = ex_lda * (1 + enhancement_correction)

    dfx_ds = 2 * mu * s / denominator**2
    exchange_gradient = ex_lda * enhancement_correction
    vx_gradient = (
        exchange_gradient
        + ex_lda / 3 * enhancement_correction
        - 4 / 3 * ex_lda * dfx_ds * s
    )
    vx_local = 4 / 3 * ex_lda + vx_gradient
    coeff_x = xp.zeros_like(norm_grad)
    nonzero_gradient = norm_grad > 0
    coeff_x[nonzero_gradient] = (
        ex_lda[nonzero_gradient]
        * dfx_ds[nonzero_gradient]
        / (
            2
            * kf[nonzero_gradient]
            * norm_grad[nonzero_gradient]
        )
    )

    # PBE correlation uses the modified Perdew-Wang LDA correlation
    # together with the PBE gradient correction H.
    ec_lda, vc_lda = lda_c_pw_mod(rho, xp)
    beta = 0.06672455060314922
    gamma = (1 - np.log(2)) / np.pi**2
    rs = (3 / (4 * np.pi * rho))**(1 / 3)
    kf_c = (9 * np.pi / 4)**(1 / 3) / rs
    ks = xp.sqrt(4 * kf_c / np.pi)
    divt = 2 * ks * rho
    t = norm_grad / divt
    exponential = xp.exp(-ec_lda / gamma)
    a = beta / (gamma * (exponential - 1))
    t2 = t**2
    at2 = a * t2
    a2t4 = at2**2
    divsum = 1 + at2 + a2t4
    div = (1 + at2) / divsum
    nolog = 1 + beta / gamma * t2 * div
    h = gamma * xp.log(nolog)

    factor = a2t4 * (2 + at2) / divsum**2
    dh = beta * t2 / nolog * (
        -7 / 3 * div
        - factor * (a * exponential * (vc_lda - ec_lda) / beta - 7 / 3)
    )
    ec = ec_lda + h
    vc_local = vc_lda + h + dh
    coeff_c = beta / (divt * ks) * (div - factor) / nolog

    gradient_field = (coeff_x + coeff_c)[..., None] * grad
    potential = vx_local + vc_local - _divergence(gradient_field, kvec, xp)
    return (
        xp.nan_to_num(ex + ec),
        xp.nan_to_num(potential),
    )


def lda_c_pw_mod(rho:np.ndarray, xp=np) -> (np.ndarray, np.ndarray):
    """
    Modified Perdew-Wang LDA correlation used by PBE.
    """
    a = 0.0310907
    a1 = 0.2137
    b1, b2, b3, b4 = 7.5957, 3.5876, 1.6382, 0.49294
    rs = (3 / (4 * np.pi * rho))**(1 / 3)
    rs12 = xp.sqrt(rs)
    rs32 = rs * rs12
    rs2 = rs**2
    omega = 2 * a * (
        b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2
    )
    logarithm = xp.log(1 + 1 / omega)
    ec = -2 * a * (1 + a1 * rs) * logarithm
    domega = 2 * a * (
        0.5 * b1 * rs12
        + b2 * rs
        + 1.5 * b3 * rs32
        + 2 * b4 * rs2
    )
    vc = (
        -2 * a * (1 + 2 / 3 * a1 * rs) * logarithm
        - 2 / 3 * a * (1 + a1 * rs)
        * domega / (omega * (omega + 1))
    )
    return ec, vc


def _gradient(field:np.ndarray, kvec:np.ndarray, xp=np) -> np.ndarray:
    """
    Calculate a periodic real-space gradient using FFT derivatives.
    """
    reciprocal = xp.fft.fftn(field)
    gradient = xp.empty((*field.shape, 3), dtype=float)
    for dim in range(3):
        gradient[..., dim] = xp.real(
            xp.fft.ifftn(1j * kvec[..., dim] * reciprocal)
        )
    return gradient


def _divergence(field:np.ndarray, kvec:np.ndarray, xp=np) -> np.ndarray:
    """
    Calculate the periodic divergence of a Cartesian vector field.
    """
    divergence = xp.zeros(field.shape[:-1], dtype=complex)
    for dim in range(3):
        divergence += 1j * kvec[..., dim] * xp.fft.fftn(field[..., dim])
    return xp.real(xp.fft.ifftn(divergence))
