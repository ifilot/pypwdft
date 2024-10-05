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

# PURPOSE
# =======
#
# Perform an electronic structure calculation of CO and visualize its valence 
# orbitals. This script builds the isosurfaces
#

import os
import numpy as np
from pytessel import PyTessel
from pypwdft import PyPWDFT, SystemBuilder, PeriodicSystem
import pickle

def main():
    # create cubic periodic system with lattice size of 10 Bohr units
    npts = 32       # number of grid points
    sz = 10
    
    if os.path.exists('ch4.pickle'):
        with open('ch4.pickle', 'rb') as f:
            res = pickle.load(f)
    else:
        # construct CO molecule system via SystemBuilder
        s = SystemBuilder().from_name('CH4', sz=sz, npts=npts)
            
        # construct calculator object
        calculator = PyPWDFT(s)
        
        # perform self-consistent field procedure and store results in res object
        res = calculator.scf(tol=1e-5, verbose=True)
        
        # print molecular orbital energies
        print(res['orbe'])
        
        # create cached result
        with open('ch4.pickle', 'wb') as f:
            pickle.dump(res, f, pickle.HIGHEST_PROTOCOL)
    
    # generate PyTessel object
    pytessel = PyTessel()
    
    # calculate overlap matrix prior to transformation
    S = calculate_overlap_matrix(res['orbc_rs'], sz, npts)
    print(S)
    
    # calculate kinetic energies prior to transformation
    for i in range(5):
        print(calculate_kinetic_energy(res['orbc_fft'][i], sz, npts).real)
    
    # perform transformation
    for i in range(5):
        res['orbc_rs'][i] = np.sign(res['orbc_rs'][i].real)*np.abs(res['orbc_rs'][i])
        
    # calculate overlap matrix after transformation
    S = calculate_overlap_matrix(res['orbc_rs'], sz, npts)
    print(S)
    
    # calculate kinetic energies after transformation
    Ct = np.sqrt(sz**3) / npts**3
    for i in range(5):
        print(calculate_kinetic_energy(np.fft.fftn(res['orbc_rs'][i]) * Ct, sz, npts).real)
    
    for i in range(5):
        print('Building isosurfaces: %02i' % (i+1))
        scalarfield = upsample_grid(res['orbc_rs'][i], sz, npts, 4)
        #scalarfield = res['orbc_rs'][i]
        unitcell = np.identity(3) * sz
        
        # build positive real isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten().real, scalarfield.shape, unitcell.flatten(), 0.03)
        pytessel.write_ply('MO_PR_%02i.ply' % (i+1), vertices, normals, indices)
        
        # build negative real isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten().real, scalarfield.shape, unitcell.flatten(), -0.03)
        pytessel.write_ply('MO_NR_%02i.ply' % (i+1), vertices, normals, indices)

def upsample_grid(scalarfield, sz, npts, upsample=4):
    """
    Upsample the grid
    """
    ct = np.sqrt(sz**3) / npts**3
    scalarfield_fft = np.fft.fftn(scalarfield) * ct
    
    Nx, Ny, Nz = scalarfield_fft.shape
    Nx_up = Nx * upsample
    Ny_up = Nx * upsample
    Nz_up = Nx * upsample

    # shift the frequencies
    fft = np.fft.fftshift(scalarfield_fft)

    # perform padding
    fft_upsampled = np.pad(fft, [((Nz_up-Nz)//2,),
                                ((Ny_up-Ny)//2,),
                                ((Nx_up-Nx)//2,)], 'constant')

    # shift back
    ct = np.sqrt(sz**3) / np.prod([Nx_up, Ny_up, Nz_up])
    fft_hires = np.fft.ifftshift(fft_upsampled) / ct

    return np.fft.ifftn(fft_hires)

def calculate_overlap_matrix(orbc, sz, npts):
    """
    Calculate the overlap matrix in real-space
    """
    N = len(orbc)
    S = np.zeros((N,N))
    dV = (sz / npts)**3
    for i in range(N):
        for j in range(N):
            S[i,j] = (np.sum(orbc[i].conjugate() * orbc[j]) * dV).real
            
    return S

def calculate_kinetic_energy(orbc_fft, sz, npts):
    """
    Calculate the kinetic energy of a molecular orbital as represented by a
    set of plane-wave coefficients
    """
    s = PeriodicSystem(sz=sz, npts=npts)
    
    return 0.5 * np.einsum('ijk,ijk,ijk', orbc_fft.conjugate(), s.get_pw_k2(), orbc_fft)

if __name__ == '__main__':
    main()