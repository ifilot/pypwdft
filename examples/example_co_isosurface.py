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

import os, sys
import numpy as np
from pytessel import PyTessel
from scipy.interpolate import RegularGridInterpolator


# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import PyPWDFT, SystemBuilder, PeriodicSystem

def main():
    # create cubic periodic system with lattice size of 10 Bohr units
    npts = 32       # number of grid points
    sz = 10
    
    # construct CO molecule system via SystemBuilder
    s = SystemBuilder().from_name('CO', sz=sz, npts=npts)
        
    # construct calculator object
    calculator = PyPWDFT(s)
    
    # perform self-consistent field procedure and store results in res object
    res = calculator.scf(tol=1e-5, nsol=9, verbose=True)
    
    # generate PyTessel object
    pytessel = PyTessel()
    
    for i in range(2,9):
        print('Building isosurfaces: %02i' % (i+1))
        scalarfield = interpolate_grid(res['orbc_rs'][i], sz, npts, 4)
        unitcell = np.identity(3) * sz
        
        # build positive real isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.real.flatten(), scalarfield.shape, unitcell.flatten(), 0.1)
        pytessel.write_ply('MO_PR_%02i.ply' % (i+1), vertices, normals, indices)
        
        # build negative real isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.real.flatten(), scalarfield.shape, unitcell.flatten(), -0.1)
        pytessel.write_ply('MO_NR_%02i.ply' % (i+1), vertices, normals, indices)
        
        # build positive imaginary isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.imag.flatten(), scalarfield.shape, unitcell.flatten(), 0.1)
        pytessel.write_ply('MO_PI_%02i.ply' % (i+1), vertices, normals, indices)
        
        # build negative imaginary isosurface
        vertices, normals, indices = pytessel.marching_cubes(scalarfield.imag.flatten(), scalarfield.shape, unitcell.flatten(), -0.1)
        pytessel.write_ply('MO_NI_%02i.ply' % (i+1), vertices, normals, indices)

def interpolate_grid(scalarfield, sz, npts, amp=2):
    x = np.linspace(0, sz, npts)
    interp = RegularGridInterpolator((x,x,x), scalarfield, method='quintic')
    s = PeriodicSystem(sz, npts * amp)
    points = s.get_r()
    
    return interp(points)

if __name__ == '__main__':
    main()