import unittest
import sys
import os
import numpy as np

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import PeriodicSystem, PyPWDFT

class TestPeriodicUnitCell(unittest.TestCase):

    def test_pwdft_ch4(self):
        """
        Test calculation of CH4
        """
        # create cubic periodic system with lattice size of 10 A and
        # 32 grid points per cartesian direction
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
        
        # test eigenvalues
        orbe = [-6.3217816, -0.5346744, -0.3280327, -0.3280327, -0.3280291]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)
        
    def test_pwdft_co(self):
        """
        Test calculation of CO
        """
        # create cubic periodic system with lattice size of 10 A and
        # 32 grid points per cartesian direction
        s = PeriodicSystem(10, 16)
        
        # add methane molecule to system
        atompos = np.array([[5.00000000, 5.00000000, 3.70973478],
                            [5.00000000, 5.00000000, 5.96769844]])
        
        # specify atomic charges
        charges = [6, 8] # C + 4 x H
        
        # add atoms to the system
        for p,c in zip(atompos, charges):
            s.add_atom(p[0], p[1], p[2], c)
        
        # perform SCF calculation
        calculator = PyPWDFT(s)
        res = calculator.scf(tol=1e-4, nsol=9, verbose=False)
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -79.69113531, 
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-7.18623545, -5.96758764, -2.65547263, -0.43149135, -0.35865529, 
                -0.35865529, -0.26187348, -0.04492418, -0.02335079]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)

    def test_pwdft_co_translation(self):
        """
        Test calculation of CO with translated molecule; the result should be
        - within numerical approximation - the same
        """
        # create cubic periodic system with lattice size of 10 A and
        # 32 grid points per cartesian direction
        s = PeriodicSystem(10, 16)
        
        # add methane molecule to system
        atompos = np.array([[5.00000000, 5.00000000, 3.70973478],
                            [5.00000000, 5.00000000, 5.96769844]])
        
        # specify atomic charges
        charges = [6, 8] # C + 4 x H
        
        # add atoms to the system
        for p,c in zip(atompos, charges):
            s.add_atom(p[0], p[1], p[2], c)
        
        # translate molecule
        s.translate((0,0,5))

        # test new positions of atoms
        np.testing.assert_almost_equal(s.get_atom_positions()[0],
                                       np.array([atompos[0,0], atompos[0,1], atompos[0,2] + 5]))
        np.testing.assert_almost_equal(s.get_atom_positions()[1],
                                       np.array([atompos[1,0], atompos[1,1], atompos[1,2] - 5]))

        # perform SCF calculation
        calculator = PyPWDFT(s)
        res = calculator.scf(tol=1e-4, nsol=9, verbose=False)
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -79.69113531, 
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-7.18623545, -5.96758764, -2.65547263, -0.43149135, -0.35865529, 
                -0.35865529, -0.26187348, -0.04492418, -0.02335079]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)

    def test_pwdft_co_npfft(self):
        """
        Test calculation of CO using Numpy FFT
        """
        # create cubic periodic system with lattice size of 10 A and
        # 32 grid points per cartesian direction
        s = PeriodicSystem(10, 16)
        
        # add methane molecule to system
        atompos = np.array([[5.00000000, 5.00000000, 3.70973478],
                            [5.00000000, 5.00000000, 5.96769844]])
        
        # specify atomic charges
        charges = [6, 8] # C + 4 x H
        
        # add atoms to the system
        for p,c in zip(atompos, charges):
            s.add_atom(p[0], p[1], p[2], c)
        
        # perform SCF calculation
        calculator = PyPWDFT(s, fft='numpy')
        res = calculator.scf(tol=1e-4, nsol=9, verbose=False)
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -79.69113531, 
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-7.18623545, -5.96758764, -2.65547263, -0.43149135, -0.35865529, 
                -0.35865529, -0.26187348, -0.04492418, -0.02335079]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)
    
    def test_pwdft_co_scipy(self):
        """
        Test calculation of CO using Numpy FFT
        """
        # create cubic periodic system with lattice size of 10 A and
        # 32 grid points per cartesian direction
        s = PeriodicSystem(10, 16)
        
        # add methane molecule to system
        atompos = np.array([[5.00000000, 5.00000000, 3.70973478],
                            [5.00000000, 5.00000000, 5.96769844]])
        
        # specify atomic charges
        charges = [6, 8] # C + 4 x H
        
        # add atoms to the system
        for p,c in zip(atompos, charges):
            s.add_atom(p[0], p[1], p[2], c)
        
        # perform SCF calculation
        calculator = PyPWDFT(s, fft='scipy')
        res = calculator.scf(tol=1e-4, nsol=9, verbose=False)
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -79.69113531, 
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-7.18623545, -5.96758764, -2.65547263, -0.43149135, -0.35865529, 
                -0.35865529, -0.26187348, -0.04492418, -0.02335079]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)

if __name__ == '__main__':
    unittest.main()
