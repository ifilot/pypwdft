from pypwdft import PyPWDFT, PeriodicSystem, SystemBuilder
import numpy as np

def main():
    # create cubic periodic system with lattice size of 10 Bohr
    ecut = 5    # wavefunction cutoff in Hartree
    sz = 10
    # construct H2 molecule system via SystemBuilder
    s = SystemBuilder().from_name('H2', sz=sz, ecut=ecut)
        
    # construct calculator object
    calculator = PyPWDFT(s)
    
    # perform self-consistent field procedure and store results in res object
    res = calculator.scf(tol=1e-6, density_tol=1e-6, verbose=True)

if __name__ == '__main__':
    main()
