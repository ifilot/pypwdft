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
# Perform an electronic structure calculation of He
#

import os, sys

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import PyPWDFT, PeriodicSystem
import numpy as np

def main():
    # create cubic periodic system with lattice size of 10 Bohr units
    npts = 16       # number of grid points
    sz = 10
    s = PeriodicSystem(sz, npts)
    
    # add helium atom to system
    atompos = np.array([[5.00000000, 5.00000000, 5.00000000]])
                        
    # specify atomic charges
    charges = [2] # He
        
    # add atoms to the system
    for p,c in zip(atompos, charges):
        s.add_atom(p[0], p[1], p[2], c)
                
    # construct calculator object
    calculator = PyPWDFT(s)
    
    # perform self-consistent field procedure and store results in res object
    res = calculator.scf(tol=1e-4, nsol=2, verbose=True)
    
    print(res['Etot'])
    print(res['Erep'])
    print(res['Enuc'])
    
if __name__ == '__main__':
    main()