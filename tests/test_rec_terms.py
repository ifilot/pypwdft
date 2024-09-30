import unittest
import sys
import os
import numpy as np

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import PeriodicSystem, PyPWDFT

class TestRecTerms(unittest.TestCase):
    """
    Test evaluation of terms in reciprocal space
    """

    def test_hartree_interaction(self):
        """
        Test equivalence assessment real and reciprocal space Hartree energy
        """
        # create cubic periodic system with lattice size of 10 A and
        # 16 grid points per cartesian direction
        s = PeriodicSystem(10, 16)
        
        # add methane molecule to system
        atompos = np.array([[5.00000000, 5.00000000, 5.00000000],
                            [6.19575624, 6.19575624, 6.19575624],
                            [3.80424376, 3.80424376, 6.19575624],
                            [3.80424376, 6.19575624, 3.80424376],
                            [6.19575624, 3.80424376, 3.80424376]])
        
        # specify atomic charges
        charges = [6, 1, 1, 1, 1] # C + 4 x H
        
        # add atoms to the system
        for p,c in zip(atompos, charges):
            s.add_atom(p[0], p[1], p[2], c)
        
        # perform SCF calculation
        calculator = PyPWDFT(s)
        res = calculator.scf(tol=1e-5, verbose=False)
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -31.60688942881068, 
                                       decimal=4)
        
        # grab total electron density, G2-scalars and cell volumes
        edens = res['edens']
        k2 = s.get_pw_k2()
        dV = res['dV']
        Omega = s.get_omega()
        npts = s.get_npts()
        
        # calculate reciprocal space charge density
        fft_edens = np.fft.fftn(edens)
        
        # solve Poisson equation in reciprocal space
        with np.errstate(divide='ignore', invalid='ignore'):
           potg = 4 * np.pi * fft_edens / k2    # reciprocal space potential
           potg[~np.isfinite(potg)] = 0         # set non-finite values to zero
           
        # convert back to real space
        harpot = np.fft.ifftn(potg)
        
        # calculate electron-electron interaction energy both in real and
        # reciprocal space
        #
        # !NOTE!: 
        # * For real fields, the FFT coefficients are complex conjugates
        #   upon G-vector inversion
        # * The Ct constant is required to obtain normalized plane waves from
        #   the FFT procedure. This constant is used twice as two sets of FFT
        #   constants are multiplied in reciprocal space.
        #
        Ct = (npts**3 / np.sqrt(Omega))
        Eee_real = np.sum(harpot * edens) * dV
        Eee = np.sum(fft_edens.conjugate() * potg) / Ct**2
        
        # assert compliance
        np.testing.assert_almost_equal(Eee_real, Eee)

if __name__ == '__main__':
    unittest.main()
