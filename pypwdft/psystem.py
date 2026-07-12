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
import math
from scipy.fft import next_fast_len

class PeriodicSystem:
    """
    Class that encapsulates a cubic unitcell with periodic boundary conditions
    and which can host the nuclei and electrons
    """
    def __init__(self, sz:float, *, ecut:float,
                 density_ecut:float=None):
        """Build a periodic system

        Args:
            sz (float): edge size of the cubic unit cell
            ecut (float): spherical wavefunction plane-wave cutoff in Hartree.
            density_ecut (float, optional): density cutoff in Hartree. Defaults
                to four times ``ecut``, as required to represent products of
                wavefunctions without aliasing.
        """
        if not np.isfinite(sz) or sz <= 0:
            raise ValueError('sz must be positive and finite.')
        if not np.isfinite(ecut) or ecut <= 0:
            raise ValueError('ecut must be positive and finite.')
        if density_ecut is None:
            density_ecut = 4 * ecut
        if not np.isfinite(density_ecut) or density_ecut < 4 * ecut:
            raise ValueError('density_ecut must be at least four times ecut.')

        # Retain both derived sizes so the solver can later use separate
        # wavefunction and density FFT grids without changing the input model.
        self.__wavefunction_npts = self.__grid_size(sz, ecut)
        self.__density_npts = self.__grid_size(sz, density_ecut)
        
        # size of the cube edges
        self.__sz = sz          
        
        # The current solver uses the density grid for all FFT operations.
        self.__npts = self.__density_npts
        
        # unit cell volume
        self.__Omega = sz**3
        
        # unit cell matrix
        self.__unitcell = np.eye(3,3) * sz
        
        # build FFT vectors and store these in the class
        self.__build_fft_vectors()

        # Select the orbital plane-wave basis.  The FFT grid is also used for
        # densities and local potentials, but only coefficients inside this
        # spherical cutoff are variational degrees of freedom.
        self.__ecut = float(ecut)
        self.__density_ecut = float(density_ecut)
        self.__pw_mask = 0.5 * self.__k2 <= self.__ecut
        
        # create placeholders for atom positions and charges
        self.__atompos = np.zeros((0,3), dtype=np.float64)
        self.__atomchg = np.array([], dtype=np.uint8)
    
    def __str__(self) -> str:
        """string representation of periodic system

        Returns:
            str: string representation
        """
        res = str(self.__unitcell) + '\n'
        for z,a in zip(self.__atomchg, self.__atompos):
            res += '%i  (%6.4f  %6.4f  %6.4f)\n' % (z,a[0], a[1], a[2])
        return res
    
    def add_atom(self, x:float, y:float, z:float, charge:float, unit:str='bohr'):
        """Add an atom to the unit cell

        Args:
            x (float): x-coordinate
            y (float): y-coordinate
            z (float): z-coordinate
            charge (float): charge of the atom
            unit (str, optional): Which length unit, 'bohr' or 'angstrom'. Defaults to 'bohr'.

        Raises:
            Exception: Invalid unit received.
        """
        if unit == 'bohr':
            pos = np.array([x,y,z], dtype=np.float64)
        elif unit == 'angstrom':
            pos = np.array([x,y,z], dtype=np.float64) * 1.8897259886 # angstrom to bohr conversion
        else:
            raise Exception('Unknown unit: %s' % unit)

        # place the atom positions at the end of the Nx3 matrix
        self.__atompos = np.vstack([self.__atompos, pos])
        
        # append the atomic charge
        self.__atomchg = np.append(self.__atomchg, charge)
    
    def get_atom_positions(self) -> np.ndarray:
        """Get atomic positions in real-space

        Returns:
            np.ndarray: atomic positions
        """
        return self.__atompos
    
    def get_atom_charges(self) -> np.ndarray:
        """Get the atomic charges

        Returns:
            np.ndarray: atomic charges
        """
        return self.__atomchg
    
    def get_omega(self) -> float:
        """Get unitcell volume

        Returns:
            float: volume of the unit cell
        """
        return self.__Omega
    
    def get_r(self) -> np.ndarray:
        """Get the sampling vectors in real-space

        Returns:
            np.ndarray: real-space sampling vectors
        """
        return self.__cvec
    
    def get_r_norms(self) -> np.ndarray:
        """Get real-space sampling vector lengths

        Returns:
            np.ndarray: real-space sampling vector lengths
        """
        return np.sqrt(np.einsum('ijkl,ijkl->ijk', self.__cvec, self.__cvec))
    
    def get_pw_k(self) -> np.ndarray:
        """Get the plane wave vectors

        Returns:
            np.ndarray: plane wave vectors
        """
        return self.__kvec
    
    def get_ct(self) -> float:
        """Get the FFT transformation constant from canonical FFT to an FFT using
           a normalized plane wave basis set.

        Returns:
            float: FFT transformation constant
        """
        return np.sqrt(self.__Omega) / self.__npts**3
    
    def get_pw_k2(self) -> np.ndarray:
        """Get the squared length of the plane wave vectors

        Returns:
            np.ndarray: squared length of plane wave vectors
        """
        return self.__k2

    def get_pw_mask(self) -> np.ndarray:
        """Get the mask selecting orbital plane waves inside ``ecut``."""
        return self.__pw_mask

    def get_ecut(self):
        """Get the wavefunction cutoff in Hartree."""
        return self.__ecut

    def get_density_ecut(self) -> float:
        """Get the density cutoff in Hartree."""
        return self.__density_ecut

    def get_wavefunction_npts(self) -> int:
        """Get the prospective standalone wavefunction FFT grid size."""
        return self.__wavefunction_npts

    def get_density_npts(self) -> int:
        """Get the density FFT grid size used by the current solver."""
        return self.__density_npts

    def get_n_plane_waves(self) -> int:
        """Get the number of plane waves in the orbital basis."""
        return int(np.count_nonzero(self.__pw_mask))
    
    def get_npts(self) -> int:
        """Get the working density-grid size per Cartesian direction.

        Returns:
            int: number of density-grid points per Cartesian direction
        """
        return self.__density_npts
    
    def get_nelec(self) -> int:
        """Get the number of electrons

        Returns:
            int: number of electrons
        """
        return np.sum(self.__atomchg)
    
    def translate(self, dist:np.ndarray):
        """Translate all atoms in unit cell

        Args:
            dist (np.ndarray): displacement vector
        """
        self.__atompos = np.mod(self.__atompos + dist, self.__sz)

    def __build_fft_vectors(self):
        """
        Construct the reciprocal space vectors of the plane waves
        """
        # determine grid points in real space
        c = np.linspace(0, self.__sz, self.__npts, endpoint=False)

        # construct real space sampling vectors
        z, y, x = np.meshgrid(c, c, c, indexing='ij')
        
        N = len(c)
        cvec = np.zeros((self.__npts,self.__npts,self.__npts,3))
        cvec[:,:,:,0] = x
        cvec[:,:,:,1] = y
        cvec[:,:,:,2] = z
        
        # calculate plane wave vector coefficients in one dimension
        k = np.fft.fftfreq(self.__npts) * 2.0 * np.pi * (self.__npts / self.__sz)
        
        # construct plane wave vectors
        k3, k2, k1 = np.meshgrid(k, k, k, indexing='ij')
        
        N = len(k)
        kvec = np.zeros((N,N,N,3))
        kvec[:,:,:,0] = k1
        kvec[:,:,:,1] = k2
        kvec[:,:,:,2] = k3
        
        k2 = np.einsum('ijkl,ijkl->ijk', kvec, kvec)
        
        self.__cvec = cvec
        self.__kvec = kvec
        self.__k2 = k2

    @staticmethod
    def __grid_size(sz:float, cutoff:float) -> int:
        """Return an FFT-friendly grid that represents ``cutoff`` exactly.

        A symmetric reciprocal grid needs indices from ``-m`` through ``m``,
        where ``m`` is the largest integer shell inside the requested cutoff.
        ``next_fast_len`` may enlarge that minimum for efficient FFTs.
        """
        gmax = np.sqrt(2 * cutoff)
        max_index = int(np.floor(sz * gmax / (2 * np.pi) + 1e-12))
        return int(next_fast_len(2 * max_index + 1))
        
    def calculate_ewald_sum(self, gcut:float=2, gamma:float=1e-8,
                            charges:np.ndarray=None) -> float:
        """Calculate Ewald sum

        Args:
            gcut (float, optional): Plane wave cut off energy in Ht. Defaults to 2.
            gamma (float, optional): Separation parameter. Defaults to 1e-8.
            charges (np.ndarray, optional): Ionic charges. Defaults to the
                stored nuclear charges.

        Returns:
            float: Ewald sum in Ht
        """
        if charges is None:
            charges = self.__atomchg
        charges = np.asarray(charges, dtype=float)
        if charges.shape != self.__atomchg.shape:
            raise ValueError('Ionic charges must match the number of atoms.')

        # establish alpha value for screening Gaussian charges
        alpha = -0.25 * gcut**2 / np.log(gamma)

        # subtract spurious self-interaction
        Eself = np.sqrt(alpha / np.pi) * np.sum(charges**2)
        
        # subtract the electroneutrality term using a uniform background charge
        Een = np.pi * np.sum(charges)**2 / (2 * alpha * self.__Omega)

        # calculate short-range interaction
        Esr = 0
        amag = np.linalg.norm(self.__unitcell, axis=1) # determine unitcell vector magnitudes
        Nmax = np.rint(np.sqrt(-0.5 * np.log(gamma)) / np.sqrt(alpha) / amag + 1.5)
        T = self.__build_indexed_vectors_excluding_zero(Nmax) @ self.__unitcell

        for ia in range(len(self.__atompos)):
            for ja in range(len(self.__atompos)):
                Rij = self.__atompos[ia] - self.__atompos[ja]       # interatomic distance
                ZiZj = charges[ia] * charges[ja]                    # product of charges
                for t in T:   # loop over all unit cell permutations
                    R = np.linalg.norm(Rij + t)
                    Esr += 0.5 * ZiZj * math.erfc(R * np.sqrt(alpha)) / R
                if ia != ja:  # terms in primary unit cell
                    R = np.linalg.norm(Rij)
                    Esr += 0.5 * ZiZj * math.erfc(R * np.sqrt(alpha)) / R

        # calculate long-range interaction
        Elr = 0
        B = 2 * np.pi * np.linalg.inv(self.__unitcell.T)            # reciprocal lattice vectors
        bm = np.linalg.norm(B, axis=1)                              # vector magnitudes
        s = np.rint(gcut / bm + 1.5)
        G = self.__build_indexed_vectors_excluding_zero(s) @ B      # produce G-vectors
        G2 = np.linalg.norm(G, axis=1)**2                           # calculate G-lengths
        pre = 2 * np.pi / self.__Omega * np.exp(-0.25 * G2 / alpha) / G2

        for ia in range(len(self.__atompos)):
            for ja in range(len(self.__atompos)):
                Rij = self.__atompos[ia] - self.__atompos[ja]
                ZiZj = charges[ia] * charges[ja]
                GR = np.sum(G * Rij, axis=1)
                Elr += ZiZj * np.sum(pre * np.cos(GR)) # discard imaginary values by using cos
        
        Eewald = Elr + Esr - Eself - Een

        return Eewald
    
    def __build_indexed_vectors_excluding_zero(self, s):
        """
        Build a set of incrementing vectors from [-s_i,s_i], exclusing the zero-term
        """
        m1 = np.arange(-s[0], s[0] + 1)
        m2 = np.arange(-s[1], s[1] + 1)
        m3 = np.arange(-s[2], s[2] + 1)
        M = np.transpose(np.meshgrid(m1, m2, m3)).reshape(-1, 3)
        return M[~np.all(M == 0, axis=1)] # remove zero-term
    
    def calculate_vpot(self) -> np.ndarray:
        """Construct the nuclear attraction potential

        Returns:
            np.ndarray: nuclear attraction potential in real-space
        """
        # calculate structure factor
        sf = np.exp(1j * self.__kvec @ self.__atompos.T)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            nucpotg = -4.0 * np.pi / self.__k2
            nucpotg[0,0,0] = 0

        # produce the nuclear attraction field      
        vpot = np.fft.fftn(np.einsum('ijk,ijkl,l->ijk', 
                                     nucpotg, 
                                     sf, 
                                     self.__atomchg)) / self.__Omega
        
        return vpot
