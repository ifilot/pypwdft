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
    def __init__(self, sys:PeriodicSystem, fft:str='pyfftw'):
        """Build PyPWDFT class

        Args:
            sys (PeriodicSystem): periodic system
            fft (str, optional): which FFT algorithm to use. Defaults to 'pyfftw'.
        """
        self.__s = sys
        self.__fft = fft
        
    def scf(self, tol:float=1e-5, nsol:int=None, verbose:bool=False) -> dict:
        """Perform self-consistent field procedure

        Args:
            tol (float, optional): electronic convergence criterion. Defaults to 1e-5.
            nsol (int, optional): number of solutions to find. Defaults to None.
            verbose (bool, optional): whether verbose output should be given. Defaults to False.

        Returns:
            dict: Dictionary containing system results
        """
        # grab a set of auxiliary variables from the PeriodicSystem class
        tstart = timeit.default_timer()     # keep track of overall time
        npts = self.__s.get_npts()          # #gridpoints per Cart. dir.
        nelec = self.__s.get_nelec()        # number of electrons
        nocc = nelec // 2                   # number of occupied orbitals
        k2 = self.__s.get_pw_k2()           # PW vector lengths
        Omega = self.__s.get_omega()        # unit cell size
        Ct = Ct = np.sqrt(Omega) / npts**3  # transformation constant
        dV = dV = Omega / npts**3           # integration constant in real space

        # number of solutions to find
        if nsol == None:       
            nsol = nocc
        else:
            if nsol < nocc:
                nsol = nocc
        
        # contruct initial search vector, this vector is kept consistent for
        # reproduction purposes
        v0 = np.eye(npts**3,1)
        
        # container for the Kohn-Sham states in reciprocal space
        mo_fft = np.zeros((nsol, npts, npts, npts), dtype=np.complex128)
        
        # containers for the Kohn-Sham states in real space
        mos = np.zeros((nsol, npts, npts, npts), dtype=np.complex128)
        
        # difference in energy between two consecutive iterative steps
        diff = 1e8
        
        # total electronic energy
        Etot = 0
        
        # iteration counter
        it = 0
        
        # pre-calculate the repulsion between the periodic set of nuclei
        Eewald = self.__s.calculate_ewald_sum()
        
        # calculate the external potential by the nuclei
        vpot = self.__s.calculate_vpot()
        
        # construct initial electron density, this density is basically
        # a homogeneous electron density that adds up to the total number
        # electrons in the unit cell
        edens = self.__build_initial_edens()
        
        # calculate initial Hartree potential
        harpot = self.__calculate_hartree_potential(edens, k2)
        
        # calculate exchange-correlation energy functional and potential
        fx, vx = self.__lda_x(edens)
        fc, vc = self.__lda_c_vwn(edens)
        fxc = fx + fc
        vxc = vx + vc
              
        # loop until electronic convergence is reached
        while diff > tol:
            # increment iteration counter
            it += 1         
            
            # store old total electronic energy
            Etotold = Etot  

            # keep track of time            
            start = timeit.default_timer()
            
            # calculate total potential
            vtot = np.real(vpot + harpot + vxc)
            
            # construct a Linear Operation class to calculate a Hamiltonian
            # element
            A = LinOpH(vtot, npts, k2, fft=self.__fft)
            e,v = scipy.sparse.linalg.eigsh(A, nsol, which='SA', v0=v0)
            
            # store new eigenvectors (Kohn-Sham states)
            for i in range(nsol):
                mo_fft[i,:,:,:] = v[:,i].reshape((npts,npts,npts))
                mos[i] = np.fft.ifftn(mo_fft[i,:,:,:]) / Ct
        
            # set mixing factor to slowly introduce the new density to the
            # old electron density
            alpha = 0.3
            
            # calculate new electron density
            edens = (1.0 - alpha) * edens + \
                    alpha * np.einsum('ijkl,ijkl->jkl', 
                                      mos[:nocc].conj(), 
                                      mos[:nocc]) * 2.0

            # calculate new Hartree potential
            harpot = self.__calculate_hartree_potential(edens, k2)
            
            # calculate new exchange potential
            fx, vx = self.__lda_x(edens)
            fc, vc = self.__lda_c_vwn(edens)
            fxc = fx + fc
            vxc = vx + vc
            
            # calculate new energy
            
            # calculate repulsion between the electrons
            Erep = np.real(0.5 * np.einsum('ijk,ijk', harpot, edens) * dV)
            
            # calculate attraction between nuclei and electrons
            Enuc = np.real(np.einsum('ijk,ijk', vpot, edens)) * dV
            
            # calculate kinetic energy (done in reciprocal space)
            Ekin = np.real(np.einsum('ijkl,ijkl,jkl', 
                           mo_fft[:nocc],
                           mo_fft[:nocc].conjugate(),
                           k2))
            
            # calculate exchange-correlation energy
            Exc = np.real(np.einsum('ijk,ijk', fxc, edens)) * dV
            
            # sum all terms to find total electronic energy
            Etot = Ekin + Enuc + Erep + Eewald + Exc
            
            # calculate difference between old and new density
            diff = np.abs(Etot - Etotold)
            
            # capture total calculation time
            stop = timeit.default_timer()
            
            # calculate time difference
            dt = stop - start
            
            # output results iteration update to user
            if verbose:
                print('%03i | Etot = %12.8f Ht | eps = %6.4e | dt = %6.4f s' \
                      % (it, Etot, diff, dt))
            
            # store lowest energy Kohn-Sham state for next iteratinon of the
            # algorithm
            v0 = v[:,0]
            
        # sort by eigenvalue in ascending order
        indx = e.argsort()
        mo_fft = mo_fft[indx]
        mos = mos[indx]
        e = e[indx]
        
        # determine total computation time
        tstop = timeit.default_timer()
        ttime = tstop - tstart
            
        # construct a dictionary that contains all the relevant output objects
        res = {
            'energy': Etot,
            'Etot': Etot,
            'Ekin': Ekin,
            'Enuc': Enuc,
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
    
    def __build_initial_edens(self):
        """
        Construct a uniform background charge to kick-start the SCF procedure
        """
        npts = self.__s.get_npts()      # number of grid points per direction
        Omega = self.__s.get_omega()    # unit cell volume
        nelec = self.__s.get_nelec()    # number of electrons
        
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
        A = 0.0621814
        x0 = -0.409286
        b = 13.0720
        c = 42.7198
    
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
    
class LinOpH(LinearOperator):
    """
    This class encapsulates the Linear Operator that functionally applies
    the Hamiltonian matrix to an input vector
    """
    def __init__(self, nu_pot:np.ndarray, npts:int, k2:np.ndarray, fft:str='pyfftw'):
        """Build Linear Operator class to solve matrix-vector multiplication in Arnoldi method

        Args:
            nu_pot (np.ndarray): real-space representation of effective potential
            npts (int): number of sampling points per Cartesian direction
            k2 (np.ndarray): squared plane wave vector lengths
            fft (str, optional): which FFT algorithm to use. Defaults to 'pyfftw'.

        Raises:
            Exception: when unknown FFT algorithm is being requested.
        """
        self.nu_pot = nu_pot
        self.npts = int(npts)
        self.k2 = k2.flatten()
        self.shape = tuple([self.npts**3, self.npts**3])
        self.dtype = np.dtype(np.complex128)
        
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
        # reshape the input (assumed to be flattened)
        psi = v.reshape((self.npts, self.npts, self.npts))
        
        # kinetic term
        T = 0.5 * self.k2 * v
        
        # return solution vector
        return T + self.__fft(psi).flatten()
    
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