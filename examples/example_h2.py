from pypwdft import PyPWDFT, PeriodicSystem, SystemBuilder
import numpy as np

def main():
    # create cubic periodic system with lattice size of 10 Bohr
    npts = 16   # number of grid points
    sz = 10
    # construct CH4 molecule system via SystemBuilder
    s = SystemBuilder().from_name('H2', sz=sz, npts=npts)
        
    # construct calculator object
    calculator = PyPWDFT(s)
    
    # perform self-consistent field procedure and store results in res object
    res = calculator.scf(tol=1e-6, density_tol=1e-6, verbose=True)

if __name__ == '__main__':
    main()