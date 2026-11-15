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
import pyfftw
import timeit
import scipy.fft
import scipy.sparse.linalg
from . import xc
from .allelectron import AllElectron
from .davidson import davidson, initial_guess
from .mixing import DensityMixer
from .psystem import PeriodicSystem

# Number of threads used by the pyFFTW plans.
FFTW_THREADS = 1

# Ratio between the eigensolver tolerance and the SCF density residual.
EIG_TOL_FACTOR = 1e-2

# Residual-norm threshold of the Davidson solver for the final, exact solve.
DAVIDSON_TOL = 1e-10


class PyPWDFT:
    """
    Class that encapsulates the planewave DFT method
    """
    def __init__(self, sys:PeriodicSystem, fft:str='pyfftw',
                 pseudopotential=None, functional:str='lda'):
        """
        Build PyPWDFT class

        Args:
            sys (PeriodicSystem): periodic system
            fft (str, optional): Numerical backend: ``'pyfftw'``, ``'numpy'``,
                ``'scipy'``, or the optional GPU backend ``'cupy'``. Defaults
                to ``'pyfftw'``.
            pseudopotential (optional): Ionic pseudopotential object. Defaults
                to an all-electron Coulomb potential.
            functional (str, optional): Exchange-correlation functional.
                Supported values are ``'lda'`` (SVWN5) and ``'pbe'``.
        """
        self.__s = sys
        self.__fft = fft.lower()
        if self.__fft not in {'pyfftw', 'numpy', 'scipy', 'cupy'}:
            raise ValueError(f'Unknown FFT backend: {fft}')
        if pseudopotential is None:
            pseudopotential = AllElectron(sys)
        elif pseudopotential.system is not sys:
            raise ValueError('The pseudopotential belongs to another system.')
        self.__ions = pseudopotential
        self.__xp = np
        self.__functional = xc.normalise(functional)

    def scf(self, tol:float=1e-5, nsol:int=None, verbose:bool=False,
            density_tol:float=None, maxiter:int=100, mixing:str='pulay',
            mixing_fraction:float=0.5, mixing_history:int=6,
            eigensolver:str='davidson') -> dict:
        """
        Perform self-consistent field procedure

        Args:
            tol (float, optional): electronic convergence criterion. Defaults to 1e-5.
            nsol (int, optional): number of solutions to find. Defaults to None.
            verbose (bool, optional): whether verbose output should be given. Defaults to False.
            density_tol (float, optional): RMS density-residual convergence
                criterion. Defaults to ``tol``.
            maxiter (int, optional): maximum number of SCF iterations.
                Defaults to 100.
            mixing (str, optional): density-mixing method, ``'pulay'`` or
                ``'linear'``. Defaults to ``'pulay'``.
            mixing_fraction (float, optional): damping applied to the density
                residual. Defaults to 0.5.
            mixing_history (int, optional): number of recent densities used by
                Pulay mixing. Defaults to 6.
            eigensolver (str, optional): ``'davidson'`` (block Davidson,
                restarted from the previous orbitals) or ``'arpack'``
                (Lanczos). Defaults to ``'davidson'``.

        Returns:
            dict: Dictionary containing system results
        """
        # grab a set of auxiliary variables from the PeriodicSystem class
        tstart = timeit.default_timer()     # keep track of overall time
        npts = self.__s.get_npts()          # #gridpoints per Cart. dir.
        nelec = self.__ions.nelec           # number of electrons
        if nelec != int(nelec):
            raise ValueError('The number of electrons must be an integer.')
        nelec = int(nelec)
        if nelec % 2:
            raise ValueError(
                'Odd-electron systems are not supported by the current '
                'closed-shell implementation.'
            )
        if nelec < 2:
            raise ValueError('At least two electrons are required.')
        if density_tol is None:
            density_tol = tol
        if tol <= 0 or density_tol <= 0:
            raise ValueError('Convergence tolerances must be positive.')
        if maxiter < 1:
            raise ValueError('maxiter must be at least one.')
        if eigensolver not in {'davidson', 'arpack'}:
            raise ValueError("eigensolver must be 'davidson' or 'arpack'.")
        nocc = nelec // 2                   # number of occupied orbitals
        gpu = self.__fft == 'cupy'
        if gpu:
            xp, sparse_linalg = _load_cupy()
        else:
            xp, sparse_linalg = np, scipy.sparse.linalg
        self.__xp = xp

        mixer = DensityMixer(
            xp,
            method=mixing,
            fraction=mixing_fraction,
            history=mixing_history,
            target_mean=nelec / self.__s.get_omega(),
        )

        k2_host = self.__s.get_pw_k2()      # PW vector lengths
        pw_mask_host = self.__s.get_pw_mask()
        k2 = xp.asarray(k2_host)
        pw_mask = xp.asarray(pw_mask_host)
        kvec = xp.asarray(self.__s.get_pw_k())
        npw = self.__s.get_n_plane_waves()  # size of the orbital basis
        Omega = self.__s.get_omega()        # unit cell size
        Ct = np.sqrt(Omega) / npts**3       # transformation constant
        dV = Omega / npts**3                # integration constant in real space

        # number of solutions to find
        if nsol == None:       
            nsol = nocc
        else:
            if nsol < nocc:
                nsol = nocc
        
        if nsol >= npw:
            raise ValueError(
                f'nsol ({nsol}) must be smaller than the number of plane '
                f'waves ({npw}).'
            )

        if verbose:
            self.__print_calculation_summary(
                nelec=nelec,
                nocc=nocc,
                nsol=nsol,
                npw=npw,
                tol=tol,
                density_tol=density_tol,
                maxiter=maxiter,
                mixing=mixer.method,
                mixing_fraction=mixer.fraction,
                mixing_history=mixer.history,
                eigensolver=eigensolver,
            )

        # construct initial search vector, this vector is kept consistent for
        # reproduction purposes
        kinetic = 0.5 * k2[pw_mask]         # kinetic energy per plane wave
        if eigensolver == 'davidson':
            v0 = initial_guess(kinetic, nsol, xp)
        elif gpu:
            v0 = xp.eye(npw, 1, dtype=xp.complex128).ravel()
        else:
            v0 = np.eye(npw, 1).ravel()
        
        # container for the Kohn-Sham states in reciprocal space
        mo_fft = xp.zeros(
            (nsol, npts, npts, npts), dtype=xp.complex128
        )
        
        # containers for the Kohn-Sham states in real space
        mos = xp.zeros((nsol, npts, npts, npts), dtype=xp.complex128)
        
        # differences between consecutive SCF steps
        diff = np.inf
        density_residual = np.inf
        
        # total electronic energy
        Etot = 0
        
        # iteration counter
        it = 0
        
        # pre-calculate the repulsion between the periodic set of ions
        Eewald = self.__ions.ionic_energy()
        
        # calculate the external potential by the ions
        vpot = xp.asarray(self.__ions.local_potential())
        nonlocal_operator = self.__ions.nonlocal_operator(xp)
        
        # construct initial electron density, this density is basically
        # a homogeneous electron density that adds up to the total number
        # electrons in the unit cell
        edens = self.__build_initial_edens(nelec)

        # The Hamiltonian operator (and its FFT plans) is built once; only its
        # effective potential changes between SCF iterations.
        hamiltonian = LinOpH(
            None, npts, k2, fft=self.__fft, pw_mask=pw_mask,
            nonlocal_operator=nonlocal_operator, xp=xp
        )
        A = sparse_linalg.LinearOperator(
            hamiltonian.shape, matvec=hamiltonian.matvec,
            dtype=hamiltonian.dtype
        )

        def solve(edens, v0, eig_tol=0):
            """
            Diagonalise the Hamiltonian of the input density ``edens`` and
            return the eigenvalues, eigenvectors, output density and the
            energy terms of that output density. ``v0`` is the eigensolver
            starting point and ``eig_tol`` its tolerance; zero requests the
            exact solve.
            """
            # Calculate the input effective potential.
            harpot = self.__calculate_hartree_potential(edens, k2)
            _, vxc = xc.evaluate(self.__functional, edens, kvec, xp)
            hamiltonian.nu_pot = xp.real(vpot + harpot + vxc)

            if eigensolver == 'davidson':
                e, v, solved = davidson(
                    hamiltonian.matmat, v0, kinetic,
                    tol=max(eig_tol, DAVIDSON_TOL), xp=xp
                )
                if eig_tol == 0 and not solved:
                    raise RuntimeError(
                        'The Davidson eigensolver did not converge.'
                    )
            else:
                e, v = sparse_linalg.eigsh(
                    A, nsol, which='SA', v0=v0, tol=eig_tol
                )

            # Sort before selecting occupied states; eigensolver ordering
            # should not be relied upon for the density construction.
            indx = e.argsort()
            e = e[indx]
            v = v[:, indx]

            # store new eigenvectors (Kohn-Sham states)
            for i in range(nsol):
                mo_fft[i].fill(0)
                mo_fft[i].reshape(-1)[pw_mask.ravel()] = v[:, i]
                mos[i] = xp.fft.ifftn(mo_fft[i]) / Ct

            output_density = xp.real(xp.einsum(
                'ijkl,ijkl->jkl', mos[:nocc].conj(), mos[:nocc]
            )) * 2.0

            # Evaluate all density-dependent energy terms using the same
            # output density as the orbitals and kinetic energy.
            harpot = self.__calculate_hartree_potential(output_density, k2)
            fxc, _ = xc.evaluate(
                self.__functional, output_density, kvec, xp
            )

            # calculate repulsion between the electrons
            Erep = xp.real(
                0.5 * xp.einsum('ijk,ijk', harpot, output_density) * dV
            )

            # calculate attraction between ions and electrons
            Enuc = xp.real(
                xp.einsum('ijk,ijk', vpot, output_density)
            ) * dV

            # calculate kinetic energy (done in reciprocal space)
            Ekin = xp.real(xp.einsum('ijkl,ijkl,jkl',
                           mo_fft[:nocc],
                           mo_fft[:nocc].conjugate(),
                           k2))

            if nonlocal_operator is None:
                Enonloc = 0.0
            else:
                applied = nonlocal_operator(v[:, :nocc])
                Enonloc = 2.0 * xp.real(
                    xp.sum(v[:, :nocc].conj() * applied)
                )

            # calculate exchange-correlation energy
            Exc = xp.real(
                xp.einsum('ijk,ijk', fxc, output_density)
            ) * dV

            energies = {
                'Ekin': self.__scalar(Ekin),
                'Enuc': self.__scalar(Enuc),
                'Enonloc': self.__scalar(Enonloc),
                'Erep': self.__scalar(Erep),
                'Exc': self.__scalar(Exc),
            }
            # sum all terms to find total electronic energy
            energies['Etot'] = sum(energies.values()) + Eewald

            return e, v, output_density, energies
        
        # loop until both energy and density are converged
        converged = False
        for it in range(1, maxiter + 1):
            # store old total electronic energy
            Etotold = Etot

            # keep track of time            
            start = timeit.default_timer()
            
            # Far from self-consistency an exact diagonalisation is wasted, so
            # the eigensolver tolerance follows the density residual. It
            # stays a factor EIG_TOL_FACTOR below both the residual and its
            # convergence threshold, which keeps eigensolver noise out of
            # the convergence test.
            eig_tol = EIG_TOL_FACTOR * min(
                max(density_residual, density_tol), 1.0
            )
            e, v, output_density, energies = solve(edens, v0, eig_tol)
            Etot = energies['Etot']
            density_residual = self.__scalar(
                xp.sqrt(xp.mean((output_density - edens)**2))
            )
            
            # calculate difference between consecutive output energies
            diff = abs(Etot - Etotold)
            
            # capture total calculation time
            stop = timeit.default_timer()
            
            # calculate time difference
            dt = stop - start
            
            # output results iteration update to user
            if verbose:
                print('%03i | Etot = %12.8f Ht | dE = %6.4e | dn = %6.4e | dt = %6.4f s' \
                      % (it, Etot, diff, density_residual, dt))
            
            # Restart the next eigensolve from the current Kohn-Sham states:
            # Davidson reuses all of them, ARPACK only the lowest one.
            v0 = v if eigensolver == 'davidson' else v[:, 0]
            if gpu and eigensolver == 'arpack':
                # CuPy's eigsh passes v0 directly to cuBLAS and mutates it.
                # A column view is strided when more than one state is
                # returned, so provide a contiguous device-owned restart.
                v0 = xp.ascontiguousarray(v0)

            if diff <= tol and density_residual <= density_tol:
                converged = True
                edens = output_density
                break

            edens = mixer.update(edens, output_density)

        if not converged:
            raise RuntimeError(
                f'SCF did not converge in {maxiter} iterations '
                f'(dE={diff:.3e}, dn={density_residual:.3e}).'
            )

        # Perform a final solve from the converged density.  This makes the
        # returned orbitals, density, and every energy component belong to one
        # consistent final output state rather than to a mixed density.
        e, v, final_density, energies = solve(edens, v0)
        density_residual = self.__scalar(
            xp.sqrt(xp.mean((final_density - edens)**2))
        )
        edens = final_density
        Etot = energies['Etot']
        
        # determine total computation time
        tstop = timeit.default_timer()
        ttime = tstop - tstart
            
        # construct a dictionary that contains all the relevant output objects
        res = {
            'energy': Etot,
            **energies,
            'edens': self.__to_numpy(edens),
            'k2': k2_host,
            'dV': dV,
            'Eewald': Eewald,
            'orbc_fft': self.__to_numpy(mo_fft), # Kohn-Sham states in reciprocal space
            'orbe': self.__to_numpy(e),          # orbital energies
            'orbc_rs': self.__to_numpy(mos),     # Kohn-Sham states in real space
            'ttime': ttime,     # total computation time
            'iterations': it,
            'energy_residual': diff,
            'density_residual': density_residual,
            'converged': converged,
            'ecut': self.__s.get_ecut(),
            'density_ecut': self.__s.get_density_ecut(),
            'wavefunction_npts': self.__s.get_wavefunction_npts(),
            'density_npts': self.__s.get_density_npts(),
            'npw': npw,
            'pseudopotential': self.__ions,
            'functional': self.__functional,
            'fft': self.__fft,
            'mixing': mixer.method,
            'mixing_fraction': mixer.fraction,
            'mixing_history': mixer.history,
            'eigensolver': eigensolver,
            'hamiltonian_applications': hamiltonian.applications,
        }
        
        return res

    def __print_calculation_summary(self, nelec, nocc, nsol, npw, tol,
                                    density_tol, maxiter, mixing,
                                    mixing_fraction, mixing_history,
                                    eigensolver):
        """
        Print the calculation setup before the SCF iterations begin.
        """
        separator = "=" * 72
        subsection = "-" * 72
        omega = self.__s.get_omega()
        cell_edge = omega ** (1.0 / 3.0)
        ecut = self.__s.get_ecut()
        density_ecut = self.__s.get_density_ecut()
        backend = {
            "cupy": "CuPy (GPU)",
            "pyfftw": "pyFFTW (CPU)",
            "numpy": "NumPy (CPU)",
            "scipy": "SciPy (CPU)",
        }[self.__fft]
        labels = self.__ions.symbols
        positions = self.__s.get_atom_positions()

        print(separator)
        print("PyPWDFT self-consistent field calculation")
        print(subsection)
        print(f"XC functional          : {self.__functional.upper()}")
        print(f"Ionic model            : {self.__ions.label}")
        print(f"FFT backend            : {backend}")
        print(
            f"Wavefunction cutoff    : {ecut:.6g} Ha "
            f"({ecut * 27.211386245988:.2f} eV)"
        )
        print(f"Density cutoff         : {density_ecut:.6g} Ha")
        print(
            f"Cubic unit cell        : {cell_edge:.6f} x {cell_edge:.6f} "
            f"x {cell_edge:.6f} bohr"
        )
        print(f"Cell volume            : {omega:.6f} bohr^3")
        print(
            f"Wavefunction FFT grid  : "
            f"{self.__s.get_wavefunction_npts()}^3"
        )
        print(
            f"Density FFT grid       : {self.__s.get_density_npts()}^3"
        )
        print(f"Plane waves            : {npw}")
        print(f"Electrons              : {nelec}")
        print(f"Occupied orbitals      : {nocc}")
        print(f"Computed orbitals      : {nsol}")
        print(f"Atoms                  : {len(positions)}")
        print(subsection)
        print("Atom positions (bohr)")
        print("  #  Atom              x             y             z")
        for index, (label, position) in enumerate(
                zip(labels, positions), start=1):
            print(
                f"{index:3d}  {label:<6s} "
                f"{position[0]:13.7f} {position[1]:13.7f} "
                f"{position[2]:13.7f}"
            )
        print(subsection)
        print(
            f"Convergence thresholds : dE <= {tol:.3e} Ha, "
            f"RMS(density) <= {density_tol:.3e}"
        )
        print(f"Maximum iterations     : {maxiter}")
        if mixing == "pulay":
            mixing_label = (
                f"Pulay (fraction={mixing_fraction:g}, "
                f"history={mixing_history})"
            )
        else:
            mixing_label = f"Linear (fraction={mixing_fraction:g})"
        print(f"Density mixing         : {mixing_label}")
        print(
            "Eigensolver            : "
            f"{'Davidson' if eigensolver == 'davidson' else 'ARPACK'}"
        )
        print(subsection)
        print("SCF iterations")
        print("Iter | Total energy (Ha) |       dE |       dn | Time (s)")
        print(subsection)

    def __scalar(self, value):
        """
        Convert a scalar value to a Python float, handling both NumPy and
        """
        if self.__xp is np:
            return float(np.real(value))
        return float(self.__xp.asnumpy(self.__xp.real(value)))

    def __to_numpy(self, value):
        """
        Convert a value to a NumPy array, handling both NumPy and CuPy arrays.
        """
        if self.__xp is np:
            return value
        return self.__xp.asnumpy(value)
    
    def __calculate_hartree_potential(self, edens:np.ndarray, k2:np.ndarray):
        """
        Calculate the Hartree potential by solving the Poisson equation
        """
        # calculate reciprocal space charge density
        xp = self.__xp
        fft_edens = xp.fft.fftn(edens)
        
        # solve Poisson equation in reciprocal space
        potg = xp.zeros_like(fft_edens)
        nonzero = k2 != 0
        potg[nonzero] = 4 * np.pi * fft_edens[nonzero] / k2[nonzero]
           
        
        # convert back to real space
        harpot = xp.fft.ifftn(potg)
        
        return harpot
    
    def __build_initial_edens(self, nelec=None):
        """
        Construct a uniform background charge to kick-start the SCF procedure
        """
        npts = self.__s.get_npts()      # number of grid points per direction
        Omega = self.__s.get_omega()    # unit cell volume
        if nelec is None:
            nelec = self.__s.get_nelec()
        
        # return a field with a homogeneous electron density that adds up
        # to the total number of electrons
        return self.__xp.ones((npts, npts, npts)) * nelec / Omega
        

class LinOpH:
    """
    Matrix-free Kohn-Sham Hamiltonian acting on plane-wave coefficients.

    The same class serves NumPy and CuPy arrays; wrap :meth:`matvec` in a
    ``LinearOperator`` of the matching library for iterative eigensolvers.
    """
    def __init__(self, nu_pot:np.ndarray, npts:int, k2:np.ndarray,
                 fft:str='pyfftw', pw_mask:np.ndarray=None,
                 nonlocal_operator=None, xp=np):
        """
        Build the Hamiltonian for matrix-vector multiplication in Arnoldi method

        Args:
            nu_pot (np.ndarray): real-space representation of effective potential
            npts (int): number of sampling points per Cartesian direction
            k2 (np.ndarray): squared plane wave vector lengths
            fft (str, optional): which FFT algorithm to use. Defaults to 'pyfftw'.
            pw_mask (np.ndarray, optional): mask selecting the spherical
                orbital plane-wave basis. Defaults to all FFT coefficients.
            nonlocal_operator (callable, optional): Function applying a
                non-local potential to active reciprocal coefficients.
            xp (module, optional): array module owning the arrays, NumPy or
                CuPy. Defaults to NumPy.

        Raises:
            Exception: when unknown FFT algorithm is being requested.
        """
        self.xp = xp
        self.nu_pot = nu_pot
        self.npts = int(npts)
        if pw_mask is None:
            pw_mask = xp.ones_like(k2, dtype=bool)
        self.pw_mask = xp.asarray(pw_mask, dtype=bool).ravel()
        self.k2 = k2.ravel()[self.pw_mask]
        self.nonlocal_operator = nonlocal_operator
        self.dtype = np.complex128
        self.shape = (len(self.k2), len(self.k2))
        self.applications = 0       # number of vectors H was applied to

        # Flat indices of the active plane waves and a persistent FFT-grid
        # buffer. Only the active entries of the buffer are ever written, so
        # the remainder stays zero and no per-call allocation is needed.
        self.__index = xp.flatnonzero(self.pw_mask)
        shape = (self.npts, self.npts, self.npts)

        if fft == 'pyfftw':
            # Planned transforms between fixed aligned buffers:
            # psi (reciprocal) -> real (real space) -> out (reciprocal).
            self.__psi = pyfftw.empty_aligned(shape, dtype='complex128')
            self.__real = pyfftw.empty_aligned(shape, dtype='complex128')
            self.__out = pyfftw.empty_aligned(shape, dtype='complex128')
            plan = dict(axes=(0, 1, 2), flags=('FFTW_MEASURE',),
                        threads=FFTW_THREADS)
            self.__ifft = pyfftw.FFTW(self.__psi, self.__real,
                                      direction='FFTW_BACKWARD', **plan)
            self.__fft = pyfftw.FFTW(self.__real, self.__out,
                                     direction='FFTW_FORWARD', **plan)
            self.__psi[:] = 0   # planning may overwrite the buffers
            self.__local = self.__local_pyfftw
        elif fft in ('numpy', 'scipy', 'cupy'):
            module = {'numpy': np.fft, 'scipy': scipy.fft, 'cupy': xp.fft}[fft]
            self.__fftn, self.__ifftn = module.fftn, module.ifftn
            self.__psi = xp.zeros(shape, dtype=np.complex128)
            self.__local = self.__local_generic
        else:
            raise Exception('Invalid FFT system: %s' % fft)

    def matvec(self, v:np.ndarray) -> np.ndarray:
        """
        Perform matrix-vector multiplication via FFTs

        Args:
            v (np.ndarray): input vector in reciprocal space

        Returns:
            np.ndarray: output result in reciprocal space
        """
        self.applications += 1

        # Embed active coefficients in the full FFT grid.
        self.__psi.reshape(-1)[self.__index] = v

        # kinetic term plus local potential term, applied in real space
        result = 0.5 * self.k2 * v + self.__local().reshape(-1)[self.__index]
        if self.nonlocal_operator is not None:
            result += self.nonlocal_operator(v)
        return result

    def matmat(self, block:np.ndarray) -> np.ndarray:
        """
        Apply the Hamiltonian to every column of a block of vectors

        Args:
            block (np.ndarray): input vectors, shape ``(npw, nstate)``

        Returns:
            np.ndarray: output vectors, shape ``(npw, nstate)``
        """
        nstate = block.shape[1]

        # The planned pyFFTW transforms act on one grid at a time.
        if self.__local == self.__local_pyfftw:
            result = np.empty_like(block)
            for i in range(nstate):
                result[:, i] = self.matvec(block[:, i])
            return result

        # Otherwise transform all states in one batched FFT, which is what
        # makes a GPU efficient.
        self.applications += nstate
        xp = self.xp
        psi = xp.zeros((nstate, self.npts**3), dtype=np.complex128)
        psi[:, self.__index] = block.T
        psi = psi.reshape((nstate, self.npts, self.npts, self.npts))
        axes = (1, 2, 3)
        local = self.__fftn(
            self.__ifftn(psi, axes=axes) * self.nu_pot, axes=axes
        )
        result = (0.5 * self.k2[:, None] * block
                  + local.reshape(nstate, -1)[:, self.__index].T)
        if self.nonlocal_operator is not None:
            result += self.nonlocal_operator(block)
        return result

    def __local_generic(self):
        """
        Apply the local potential to the buffered state; return it on the
        full reciprocal FFT grid.
        """
        return self.__fftn(self.__ifftn(self.__psi) * self.nu_pot)

    def __local_pyfftw(self):
        """
        Same as :meth:`__local_generic`, using the planned in-place buffers.
        """
        self.__ifft()
        np.multiply(self.__real, self.nu_pot, out=self.__real)
        return self.__fft()

def _load_cupy():
    """Load the optional CuPy backend and its sparse linear algebra module."""
    try:
        import cupy as cp
        import cupyx.scipy.sparse.linalg as cupy_sparse_linalg
    except ImportError as exc:
        raise ImportError(
            "fft='cupy' requires a CuPy package compatible with the installed "
            "CUDA toolkit. See https://docs.cupy.dev/en/stable/install.html."
        ) from exc
    try:
        if cp.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError("fft='cupy' requires an available CUDA GPU.")
    except cp.cuda.runtime.CUDARuntimeError as exc:
        raise RuntimeError(
            "CuPy is installed, but no usable CUDA GPU could be initialized."
        ) from exc
    return cp, cupy_sparse_linalg
