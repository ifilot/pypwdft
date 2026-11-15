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
FFT backends and thread management.

:class:`GridFFT` gives every part of the solver (Hartree potential, GGA
derivatives, orbital back-transforms) the same ``fftn``/``ifftn`` interface as
``numpy.fft``, so that on the CPU these transforms use the planned, threaded
pyFFTW path instead of NumPy's single-threaded one.
"""

import contextlib
import os
import warnings

import numpy as np
import pyfftw
import scipy.fft

try:
    import threadpoolctl
except ImportError:                                     # pragma: no cover
    threadpoolctl = None


def resolve_threads(threads) -> int:
    """
    Validate a requested FFT thread count and return it as an integer.
    """
    if isinstance(threads, bool) or int(threads) != threads or threads < 1:
        raise ValueError('threads must be a positive integer.')
    return int(threads)


def openmp_thread_limit():
    """
    Return the number of threads OpenMP will hand out, or ``None`` if unknown.
    """
    limits = []
    if threadpoolctl is not None:
        limits = [
            info['num_threads'] for info in threadpoolctl.threadpool_info()
            if info.get('user_api') == 'openmp'
        ]
    for variable in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT'):
        try:
            limits.append(int(os.environ[variable].split(',')[0]))
        except (KeyError, ValueError):
            pass
    return min(limits) if limits else None


def check_fftw_threads(threads:int) -> None:
    """
    Warn when OpenMP settings prevent pyFFTW from using ``threads`` threads.

    pyFFTW builds that thread through OpenMP silently obey ``OMP_NUM_THREADS``
    rather than the thread count requested for the plan.
    """
    if threads == 1 or getattr(pyfftw, '_threading_type', None) != 'OMP':
        return
    limit = openmp_thread_limit()
    if limit is not None and limit < threads:
        warnings.warn(
            f'{threads} FFT threads were requested, but OpenMP is limited to '
            f'{limit} thread(s) (OMP_NUM_THREADS); pyFFTW will not use more. '
            'Unset OMP_NUM_THREADS or raise it to at least the requested '
            'thread count.',
            RuntimeWarning,
            stacklevel=3,
        )


def limit_blas_threads(fft_threads:int):
    """
    Context manager restricting BLAS to one thread while FFTs are threaded.

    Threaded FFTW and a threaded BLAS compete for the same cores: the idle
    threads of one library keep spinning while the other is working, which
    was measured to slow an SCF run down several-fold. The limit has to span
    the whole SCF procedure: restricting it to the FFT calls alone was
    measured to leave most of the slowdown in place. With a single FFT
    thread the BLAS thread count is left untouched.
    """
    if fft_threads == 1 or threadpoolctl is None:
        return contextlib.nullcontext()
    return threadpoolctl.threadpool_limits(limits=1, user_api='blas')


class GridFFT:
    """
    ``fftn``/``ifftn`` over whole arrays for one numerical backend.

    Args:
        backend (str): ``'pyfftw'``, ``'numpy'``, ``'scipy'`` or ``'cupy'``.
        xp (module, optional): array module owning the data, NumPy or CuPy.
        threads (int, optional): number of threads for the pyFFTW and SciPy
            backends. Defaults to 1.
    """
    def __init__(self, backend:str, xp=np, threads:int=1):
        self.backend = backend
        self.threads = resolve_threads(threads)
        if backend == 'pyfftw':
            # Plans are created on first use and reused for every later
            # transform of the same grid shape.
            self.__plans = {}
            self.fftn = self.__pyfftw_forward
            self.ifftn = self.__pyfftw_backward
        elif backend == 'scipy':
            self.fftn = self.__scipy_forward
            self.ifftn = self.__scipy_backward
        elif backend in ('numpy', 'cupy'):
            self.fftn, self.ifftn = xp.fft.fftn, xp.fft.ifftn
        else:
            raise ValueError(f'Unknown FFT backend: {backend}')

    def __plan(self, shape, direction):
        """Return ``(input buffer, plan)`` for a grid shape and direction."""
        key = (shape, direction)
        if key not in self.__plans:
            source = pyfftw.empty_aligned(shape, dtype='complex128')
            target = pyfftw.empty_aligned(shape, dtype='complex128')
            plan = pyfftw.FFTW(
                source, target, axes=tuple(range(len(shape))),
                direction=direction, flags=('FFTW_MEASURE',),
                threads=self.threads,
            )
            self.__plans[key] = (source, plan)
        return self.__plans[key]

    def __pyfftw_forward(self, field):
        source, plan = self.__plan(np.shape(field), 'FFTW_FORWARD')
        source[...] = field
        return plan().copy()

    def __pyfftw_backward(self, field):
        source, plan = self.__plan(np.shape(field), 'FFTW_BACKWARD')
        source[...] = field
        return plan().copy()

    def __scipy_forward(self, field):
        return scipy.fft.fftn(field, workers=self.threads)

    def __scipy_backward(self, field):
        return scipy.fft.ifftn(field, workers=self.threads)
