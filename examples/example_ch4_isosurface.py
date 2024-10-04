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

import numpy as np
from pytessel import PyTessel
from pypwdft import PyPWDFT, SystemBuilder

def main():
    # create cubic periodic system with lattice size of 10 Bohr units
    npts = 32       # number of grid points
    sz = 10
    
    # construct CO molecule system via SystemBuilder
    s = SystemBuilder().from_name('CH4', sz=sz, npts=npts)
        
    # construct calculator object
    calculator = PyPWDFT(s)
    
    # perform self-consistent field procedure and store results in res object
    res = calculator.scf(tol=1e-5, verbose=True)
    
    # print molecular orbital energies
    print(res['orbe'])
    
    # generate PyTessel object
    pytessel = PyTessel()
    
    for i in range(5):
        print('Building isosurfaces: %02i' % (i+1))
        scalarfield = upsample_grid(res['orbc_fft'][i], sz**3, 4)
        unitcell = np.identity(3) * sz
        
        # build positive real isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.real.flatten(), scalarfield.shape, unitcell.flatten(), 0.03)
        pytessel.write_ply('MO_PR_%02i.ply' % (i+1), vertices, normals, indices)
        
        # build negative real isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.real.flatten(), scalarfield.shape, unitcell.flatten(), -0.03)
        pytessel.write_ply('MO_NR_%02i.ply' % (i+1), vertices, normals, indices)
        
        # build positive imaginary isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.imag.flatten(), scalarfield.shape, unitcell.flatten(), 0.03)
        pytessel.write_ply('MO_PI_%02i.ply' % (i+1), vertices, normals, indices)
        
        # build negative imaginary isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.imag.flatten(), scalarfield.shape, unitcell.flatten(), -0.03)
        pytessel.write_ply('MO_NI_%02i.ply' % (i+1), vertices, normals, indices)

def upsample_grid(scalarfield_fft, Omega, upsample=4):
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
    fft_hires = np.fft.ifftshift(fft_upsampled)

    return np.fft.ifftn(fft_hires) * np.prod([Nx_up, Ny_up, Nz_up]) / np.sqrt(Omega)

if __name__ == '__main__':
    main()