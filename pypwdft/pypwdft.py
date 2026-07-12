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
from scipy.sparse.linalg import LinearOperator
import scipy.sparse.linalg
from .psystem import PeriodicSystem

class PyPWDFT:
    """
    Class that encapsulates the planewave DFT method
    """
    def __init__(self, sys:PeriodicSystem, fft:str='pyfftw',
                 pseudopotential=None, functional:str='lda'):
        """Build PyPWDFT class

        Args:
            sys (PeriodicSystem): periodic system
            fft (str, optional): which FFT algorithm to use. Defaults to 'pyfftw'.
            pseudopotential (optional): Ionic pseudopotential object. Defaults
                to an all-electron Coulomb potential.
            functional (str, optional): Exchange-correlation functional.
                Supported values are ``'lda'`` (SVWN5) and ``'pbe'``.
        """
        self.__s = sys
        self.__fft = fft
        self.__pseudopotential = pseudopotential
        functional = functional.lower()
        aliases = {'svwn5': 'lda'}
        self.__functional = aliases.get(functional, functional)
        if self.__functional not in {'lda', 'pbe'}:
            raise ValueError(
                f'Unknown exchange-correlation functional: {functional}'
            )
        if (pseudopotential is not None and
                pseudopotential.system is not sys):
            raise ValueError('The pseudopotential belongs to another system.')
        
    def scf(self, tol:float=1e-5, nsol:int=None, verbose:bool=False,
            density_tol:float=None, maxiter:int=100) -> dict:
        """Perform self-consistent field procedure

        Args:
            tol (float, optional): electronic convergence criterion. Defaults to 1e-5.
            nsol (int, optional): number of solutions to find. Defaults to None.
            verbose (bool, optional): whether verbose output should be given. Defaults to False.
            density_tol (float, optional): RMS density-residual convergence
                criterion. Defaults to ``tol``.
            maxiter (int, optional): maximum number of SCF iterations.
                Defaults to 100.

        Returns:
            dict: Dictionary containing system results
        """
        # grab a set of auxiliary variables from the PeriodicSystem class
        tstart = timeit.default_timer()     # keep track of overall time
        npts = self.__s.get_npts()          # #gridpoints per Cart. dir.
        if self.__pseudopotential is None:
            nelec = self.__s.get_nelec()    # number of electrons
        else:
            nelec = self.__pseudopotential.nelec
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
        nocc = nelec // 2                   # number of occupied orbitals
        k2 = self.__s.get_pw_k2()           # PW vector lengths
        pw_mask = self.__s.get_pw_mask()    # selected orbital plane waves
        npw = self.__s.get_n_plane_waves()  # size of the orbital basis
        Omega = self.__s.get_omega()        # unit cell size
        Ct = Ct = np.sqrt(Omega) / npts**3  # transformation constant
        dV = dV = Omega / npts**3           # integration constant in real space

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

        # construct initial search vector, this vector is kept consistent for
        # reproduction purposes
        v0 = np.eye(npw, 1)
        
        # container for the Kohn-Sham states in reciprocal space
        mo_fft = np.zeros((nsol, npts, npts, npts), dtype=np.complex128)
        
        # containers for the Kohn-Sham states in real space
        mos = np.zeros((nsol, npts, npts, npts), dtype=np.complex128)
        
        # differences between consecutive SCF steps
        diff = np.inf
        density_residual = np.inf
        
        # total electronic energy
        Etot = 0
        
        # iteration counter
        it = 0
        
        # pre-calculate the repulsion between the periodic set of nuclei
        if self.__pseudopotential is None:
            Eewald = self.__s.calculate_ewald_sum()
        else:
            Eewald = self.__pseudopotential.ionic_energy()
        
        # calculate the external potential by the nuclei
        if self.__pseudopotential is None:
            vpot = self.__s.calculate_vpot()
            nonlocal_operator = None
        else:
            vpot = self.__pseudopotential.local_potential()
            nonlocal_operator = self.__pseudopotential.apply_nonlocal
        
        # construct initial electron density, this density is basically
        # a homogeneous electron density that adds up to the total number
        # electrons in the unit cell
        edens = self.__build_initial_edens(nelec)
        
        # loop until both energy and density are converged
        converged = False
        for it in range(1, maxiter + 1):
            # store old total electronic energy
            Etotold = Etot

            # keep track of time            
            start = timeit.default_timer()
            
            # Calculate the input effective potential.
            harpot = self.__calculate_hartree_potential(edens, k2)
            fxc, vxc = self.__calculate_xc(edens)
            vtot = np.real(vpot + harpot + vxc)
            
            # construct a Linear Operation class to calculate a Hamiltonian
            # element
            A = LinOpH(vtot, npts, k2, fft=self.__fft,
                       pw_mask=pw_mask,
                       nonlocal_operator=nonlocal_operator)
            e,v = scipy.sparse.linalg.eigsh(A, nsol, which='SA', v0=v0)

            # Sort before selecting occupied states; eigensolver ordering
            # should not be relied upon for the density construction.
            indx = e.argsort()
            e = e[indx]
            v = v[:, indx]
            
            # store new eigenvectors (Kohn-Sham states)
            for i in range(nsol):
                mo_fft[i].fill(0)
                mo_fft[i].flat[pw_mask.ravel()] = v[:,i]
                mos[i] = np.fft.ifftn(mo_fft[i,:,:,:]) / Ct

            output_density = np.real(np.einsum(
                'ijkl,ijkl->jkl', mos[:nocc].conj(), mos[:nocc]
            )) * 2.0
            density_residual = np.sqrt(np.mean((output_density - edens)**2))

            # set mixing factor to slowly introduce the new density to the
            # old electron density
            alpha = 0.3

            # Evaluate all density-dependent energy terms using the same
            # output density as the orbitals and kinetic energy.
            harpot = self.__calculate_hartree_potential(output_density, k2)
            fxc, _ = self.__calculate_xc(output_density)
            
            # calculate repulsion between the electrons
            Erep = np.real(0.5 * np.einsum('ijk,ijk', harpot, output_density) * dV)
            
            # calculate attraction between nuclei and electrons
            Enuc = np.real(np.einsum('ijk,ijk', vpot, output_density)) * dV
            
            # calculate kinetic energy (done in reciprocal space)
            Ekin = np.real(np.einsum('ijkl,ijkl,jkl',
                           mo_fft[:nocc],
                           mo_fft[:nocc].conjugate(),
                           k2))

            if self.__pseudopotential is None:
                Enonloc = 0.0
            else:
                Enonloc = self.__pseudopotential.nonlocal_energy(
                    v[:, :nocc]
                )
            
            # calculate exchange-correlation energy
            Exc = np.real(np.einsum('ijk,ijk', fxc, output_density)) * dV
            
            # sum all terms to find total electronic energy
            Etot = Ekin + Enuc + Enonloc + Erep + Eewald + Exc
            
            # calculate difference between consecutive output energies
            diff = np.abs(Etot - Etotold)
            
            # capture total calculation time
            stop = timeit.default_timer()
            
            # calculate time difference
            dt = stop - start
            
            # output results iteration update to user
            if verbose:
                print('%03i | Etot = %12.8f Ht | dE = %6.4e | dn = %6.4e | dt = %6.4f s' \
                      % (it, Etot, diff, density_residual, dt))
            
            # store lowest energy Kohn-Sham state for next iteratinon of the
            # algorithm
            v0 = v[:,0]

            if diff <= tol and density_residual <= density_tol:
                converged = True
                edens = output_density
                break

            edens = (1.0 - alpha) * edens + alpha * output_density

        if not converged:
            raise RuntimeError(
                f'SCF did not converge in {maxiter} iterations '
                f'(dE={diff:.3e}, dn={density_residual:.3e}).'
            )

        # Perform a final solve from the converged density.  This makes the
        # returned orbitals, density, and every energy component belong to one
        # consistent final output state rather than to a mixed density.
        harpot = self.__calculate_hartree_potential(edens, k2)
        _, vxc = self.__calculate_xc(edens)
        vtot = np.real(vpot + harpot + vxc)
        A = LinOpH(vtot, npts, k2, fft=self.__fft,
                   pw_mask=pw_mask,
                   nonlocal_operator=nonlocal_operator)
        e, v = scipy.sparse.linalg.eigsh(A, nsol, which='SA', v0=v0)
        indx = e.argsort()
        e = e[indx]
        v = v[:, indx]
        for i in range(nsol):
            mo_fft[i].fill(0)
            mo_fft[i].flat[pw_mask.ravel()] = v[:, i]
            mos[i] = np.fft.ifftn(mo_fft[i]) / Ct
        final_density = np.real(np.einsum(
            'ijkl,ijkl->jkl', mos[:nocc].conj(), mos[:nocc]
        )) * 2.0
        density_residual = np.sqrt(np.mean((final_density - edens)**2))
        edens = final_density

        harpot = self.__calculate_hartree_potential(edens, k2)
        fxc, _ = self.__calculate_xc(edens)
        Erep = np.real(0.5 * np.einsum('ijk,ijk', harpot, edens) * dV)
        Enuc = np.real(np.einsum('ijk,ijk', vpot, edens)) * dV
        Ekin = np.real(np.einsum(
            'ijkl,ijkl,jkl', mo_fft[:nocc], mo_fft[:nocc].conjugate(), k2
        ))
        if self.__pseudopotential is None:
            Enonloc = 0.0
        else:
            Enonloc = self.__pseudopotential.nonlocal_energy(v[:, :nocc])
        Exc = np.real(np.einsum('ijk,ijk', fxc, edens)) * dV
        Etot = Ekin + Enuc + Enonloc + Erep + Eewald + Exc
        
        # determine total computation time
        tstop = timeit.default_timer()
        ttime = tstop - tstart
            
        # construct a dictionary that contains all the relevant output objects
        res = {
            'energy': Etot,
            'Etot': Etot,
            'Ekin': Ekin,
            'Enuc': Enuc,
            'Enonloc': Enonloc,
            'Erep': Erep,
            'Exc': Exc,
            'edens': edens,
            'k2': k2,
            'dV': dV,
            'Eewald': Eewald,
            'orbc_fft': mo_fft, # Kohn-Sham states in reciprocal space
            'orbe': e,          # orbital energies
            'orbc_rs': mos,     # Kohn-Sham states in real space
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
            'pseudopotential': self.__pseudopotential,
            'functional': self.__functional,
        }
        
        return res
    
    def __calculate_hartree_potential(self, edens:np.ndarray, k2:np.ndarray):
        """
        Calculate the Hartree potential by solving the Poisson equation
        """
        # calculate reciprocal space charge density
        fft_edens = np.fft.fftn(edens)
        
        # solve Poisson equation in reciprocal space
        with np.errstate(divide='ignore', invalid='ignore'):
           potg = 4 * np.pi * fft_edens / k2    # reciprocal space potential
           potg[~np.isfinite(potg)] = 0         # set non-finite values to zero
           
        
        # convert back to real space
        harpot = np.fft.ifftn(potg)
        
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
        return np.ones((npts, npts, npts)) * nelec / Omega
        
    def __lda_x(self, rho:np.ndarray) -> (np.ndarray, np.ndarray):
        """
        Slater exchange functional, see Parr and Yang pages 154 and 155,
        equations 7.4.5 and 7.4.9
        """
        f = -3 / 4 * (3 / (2 * np.pi))**(2 / 3)
        rs = (3 / (4 * np.pi * rho))**(1 / 3)

        ex = f / rs
        vx = 4 / 3 * ex

        return ex, vx
    
    def __lda_c_vwn(self, rho:np.ndarray) -> (np.ndarray, np.ndarray):
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
        atan = np.arctan(Q / (2*x+b))
    
        ec = A/2 * (np.log(x**2/X) + 2*b/Q * atan - b*x0/X0 * \
                    (np.log((x-x0)**2 / X) + 2 * (b + 2*x0) / Q * atan))
    
        tt = tx**2 + Q**2
        vc = ec - x * A / 12 * (2 / x - tx / X - 4 * b / tt - fx0 * \
            (2 / (x - x0) - tx / X - 4 * (2 * x0 + b) / tt))
    
        return ec,vc

    def __calculate_xc(self, rho:np.ndarray) -> (np.ndarray, np.ndarray):
        """Evaluate XC energy per particle and its multiplicative potential."""
        if self.__functional == 'lda':
            ex, vx = self.__lda_x(rho)
            ec, vc = self.__lda_c_vwn(rho)
            return ex + ec, vx + vc
        return self.__pbe_xc(rho)

    def __pbe_xc(self, rho:np.ndarray) -> (np.ndarray, np.ndarray):
        """Spin-unpolarized PBE energy per particle and GGA potential."""
        rho = np.maximum(np.asarray(rho, dtype=float), 1e-14)
        grad = self.__gradient(rho)
        norm_grad = np.linalg.norm(grad, axis=-1)

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
        with np.errstate(divide='ignore', invalid='ignore'):
            coeff_x = np.where(
                norm_grad > 0,
                ex_lda * dfx_ds / (2 * kf * norm_grad),
                0,
            )

        # PBE correlation uses the modified Perdew-Wang LDA correlation
        # together with the PBE gradient correction H.
        ec_lda, vc_lda = self.__lda_c_pw_mod(rho)
        beta = 0.06672455060314922
        gamma = (1 - np.log(2)) / np.pi**2
        rs = (3 / (4 * np.pi * rho))**(1 / 3)
        kf_c = (9 * np.pi / 4)**(1 / 3) / rs
        ks = np.sqrt(4 * kf_c / np.pi)
        divt = 2 * ks * rho
        t = norm_grad / divt
        exponential = np.exp(-ec_lda / gamma)
        a = beta / (gamma * (exponential - 1))
        t2 = t**2
        at2 = a * t2
        a2t4 = at2**2
        divsum = 1 + at2 + a2t4
        div = (1 + at2) / divsum
        nolog = 1 + beta / gamma * t2 * div
        h = gamma * np.log(nolog)

        factor = a2t4 * (2 + at2) / divsum**2
        dh = beta * t2 / nolog * (
            -7 / 3 * div
            - factor * (a * exponential * (vc_lda - ec_lda) / beta - 7 / 3)
        )
        ec = ec_lda + h
        vc_local = vc_lda + h + dh
        coeff_c = beta / (divt * ks) * (div - factor) / nolog

        gradient_field = (coeff_x + coeff_c)[..., None] * grad
        potential = vx_local + vc_local - self.__divergence(gradient_field)
        return (
            np.nan_to_num(ex + ec),
            np.nan_to_num(potential),
        )

    def __lda_c_pw_mod(self, rho:np.ndarray) -> (np.ndarray, np.ndarray):
        """Modified Perdew-Wang LDA correlation used by PBE."""
        a = 0.0310907
        a1 = 0.2137
        b1, b2, b3, b4 = 7.5957, 3.5876, 1.6382, 0.49294
        rs = (3 / (4 * np.pi * rho))**(1 / 3)
        rs12 = np.sqrt(rs)
        rs32 = rs * rs12
        rs2 = rs**2
        omega = 2 * a * (
            b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2
        )
        logarithm = np.log(1 + 1 / omega)
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

    def __gradient(self, field:np.ndarray) -> np.ndarray:
        """Calculate a periodic real-space gradient using FFT derivatives."""
        reciprocal = np.fft.fftn(field)
        kvec = self.__s.get_pw_k()
        gradient = np.empty((*field.shape, 3), dtype=float)
        for dim in range(3):
            gradient[..., dim] = np.real(
                np.fft.ifftn(1j * kvec[..., dim] * reciprocal)
            )
        return gradient

    def __divergence(self, field:np.ndarray) -> np.ndarray:
        """Calculate the periodic divergence of a Cartesian vector field."""
        kvec = self.__s.get_pw_k()
        divergence = np.zeros(field.shape[:-1], dtype=complex)
        for dim in range(3):
            divergence += 1j * kvec[..., dim] * np.fft.fftn(field[..., dim])
        return np.real(np.fft.ifftn(divergence))
    
class LinOpH(LinearOperator):
    """
    This class encapsulates the Linear Operator that functionally applies
    the Hamiltonian matrix to an input vector
    """
    def __init__(self, nu_pot:np.ndarray, npts:int, k2:np.ndarray,
                 fft:str='pyfftw', pw_mask:np.ndarray=None,
                 nonlocal_operator=None):
        """Build Linear Operator class to solve matrix-vector multiplication in Arnoldi method

        Args:
            nu_pot (np.ndarray): real-space representation of effective potential
            npts (int): number of sampling points per Cartesian direction
            k2 (np.ndarray): squared plane wave vector lengths
            fft (str, optional): which FFT algorithm to use. Defaults to 'pyfftw'.
            pw_mask (np.ndarray, optional): mask selecting the spherical
                orbital plane-wave basis. Defaults to all FFT coefficients.
            nonlocal_operator (callable, optional): Function applying a
                non-local potential to active reciprocal coefficients.

        Raises:
            Exception: when unknown FFT algorithm is being requested.
        """
        self.nu_pot = nu_pot
        self.npts = int(npts)
        if pw_mask is None:
            pw_mask = np.ones_like(k2, dtype=bool)
        self.pw_mask = np.asarray(pw_mask, dtype=bool).flatten()
        self.k2 = k2.flatten()[self.pw_mask]
        self.nonlocal_operator = nonlocal_operator
        super().__init__(
            dtype=np.complex128,
            shape=(len(self.k2), len(self.k2)),
        )
        
        if fft == 'pyfftw':
            self.fft_in = pyfftw.empty_aligned(
                (self.npts, self.npts, self.npts), 
                dtype='complex128')
            self.fft_obj = pyfftw.builders.fftn(self.fft_in, 
                                                auto_align_input = True, 
                                                auto_contiguous = True,
                                                overwrite_input = False, 
                                                avoid_copy = False,
                                                planner_effort='FFTW_MEASURE')
            self.ifft_obj = pyfftw.builders.ifftn(self.fft_in,
                                                  auto_align_input = True, 
                                                  auto_contiguous = True,
                                                  overwrite_input = False, 
                                                  avoid_copy = False,
                                                  planner_effort='FFTW_MEASURE')
            self.__fft = self._solve_pyfft
            
        elif fft == 'numpy':
            self.__fft = self._solve_npfft
        elif fft == 'scipy':
            self.__fft = self._solve_scfft
        else:
            raise Exception('Invalid FFT system: %s' % fft)
    
    def _matvec(self, v:np.ndarray) -> np.ndarray:
        """Perform matrix-vector multiplication via FFTs

        Args:
            v (np.ndarray): input vector in reciprocal space

        Returns:
            np.ndarray: output result in reciprocal space
        """
        # Embed active coefficients in the full FFT grid.
        psi = np.zeros(self.npts**3, dtype=np.complex128)
        psi[self.pw_mask] = v
        psi = psi.reshape((self.npts, self.npts, self.npts))
        
        # kinetic term
        T = 0.5 * self.k2 * v
        
        # return solution vector
        result = T + self.__fft(psi).flatten()[self.pw_mask]
        if self.nonlocal_operator is not None:
            result += self.nonlocal_operator(v)
        return result
    
    def _solve_npfft(self, psi:np.ndarray) -> np.ndarray:
        """Solve potential interaction using Numpy FFT

        Args:
            psi (np.ndarray): molecular orbital coefficients in reciprocal-space

        Returns:
            np.ndarray: result in reciprocal-space
        """
        return np.fft.fftn(np.fft.ifftn(psi) * self.nu_pot)
        
    def _solve_pyfft(self, psi:np.ndarray) -> np.ndarray:
        """Solve potential interaction using PyFFT

        Args:
            psi (np.ndarray): molecular orbital coefficients in reciprocal-space

        Returns:
            np.ndarray: result in reciprocal-space
        """
        return self.fft_obj(self.ifft_obj(psi) * self.nu_pot)
    
    def _solve_scfft(self, psi:np.ndarray) -> np.ndarray:
        """Solve potential interaction using Scipy FFT

        Args:
            psi (np.ndarray): molecular orbital coefficients in reciprocal-space

        Returns:
            np.ndarray: result in reciprocal-space
        """
        return scipy.fft.fftn(scipy.fft.ifftn(psi) * self.nu_pot)
