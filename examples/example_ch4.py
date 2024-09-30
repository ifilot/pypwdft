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
# Perform an electronic structure calculation of CH4 and visualize its occupied
# orbitals. This script uses a low number of grid points and a very loose 
# convergence tolerance such that it is relatively fast to run.
#

import os, sys

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import PyPWDFT, PeriodicSystem, SystemBuilder
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():
    # create cubic periodic system with lattice size of 10 Bohr units
    npts = 16   # number of grid points
    sz = 10
    
    # construct CH4 molecule system via SystemBuilder
    s = SystemBuilder().from_name('CH4', sz=sz, npts=npts)
        
    # construct calculator object
    calculator = PyPWDFT(s)
    
    # perform self-consistent field procedure and store results in res object
    res = calculator.scf(tol=1e-1, verbose=True)
    
    # visualize the occupied molecular orbitals
    fig, im = plt.subplots(3,5, dpi=300, figsize=(16,8))
    extent=[-sz/2,sz/2,-sz/2,sz/2]
    orbe = res['orbe']
    fig.suptitle('Occupied molecular orbitals of CH$_{4}$')
    m = np.empty((3,5), dtype=object) # create placeholder for maps
    
    for i in range(0,5):
        # visualize the real part of the wave function
        field = np.real(res['orbc_rs'][i][npts//2, :, :])
        maxval = max(np.max(np.abs(field)), 0.1)
        m[0][i] = im[0,i].imshow(field, origin='lower',
                   interpolation='bicubic', extent=extent, cmap='PRGn',
                   vmin=-maxval, vmax=maxval)
        im[0,i].set_title(r'MO($\mathbb{R}$)%i / E=%6.4f Ht' % (i+1, orbe[i]))
        
        # visualize the imaginary part of the wave function
        field = np.imag(res['orbc_rs'][i][npts//2, :, :])
        maxval = max(np.max(np.abs(field)), 0.001)
        m[1][i] = im[1,i].imshow(field, origin='lower',
                   interpolation='bicubic', extent=extent, cmap='PRGn',
                   vmin=-maxval, vmax=maxval)
        im[1,i].set_title(r'MO($\mathbb{I}$)%i / E=%6.4f Ht' % (i+1, orbe[i]))
        
        # visualize the electron density
        m[2][i] = im[2,i].imshow(np.real(res['orbc_rs'][i][npts//2, :, :].conj() * 
                               res['orbc_rs'][i][npts//2, :, :]), 
                       origin='lower', interpolation='bicubic', extent=extent)
        im[2,i].set_title(r'MO($\rho$)%i / E=%6.4f Ht' % (i+1, orbe[i]))
        
    for j in range(0,3):
        for i in range(0,5):
            im[j,i].set_xlabel('$x$-position [$\AA$]')
            im[j,i].set_ylabel('$y$-position [$\AA$]')
            
            divider = make_axes_locatable(im[j,i])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(m[j][i], cax=cax, orientation='vertical')
        
    plt.tight_layout()
    
if __name__ == '__main__':
    main()