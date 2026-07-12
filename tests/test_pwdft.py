import unittest
import sys
import os
import numpy as np
import pytest

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import GTHPseudopotential, PeriodicSystem, PyPWDFT

class TestPeriodicUnitCell(unittest.TestCase):

    def test_odd_electron_system_is_rejected(self):
        s = PeriodicSystem(10, ecut=1)
        s.add_atom(5, 5, 5, 1)

        with self.assertRaisesRegex(ValueError, 'Odd-electron'):
            PyPWDFT(s, fft='numpy').scf()

    def test_scf_iteration_limit(self):
        s = PeriodicSystem(10, ecut=1)
        s.add_atom(4.3, 5, 5, 1)
        s.add_atom(5.7, 5, 5, 1)

        with self.assertRaisesRegex(RuntimeError, 'did not converge'):
            PyPWDFT(s, fft='numpy').scf(maxiter=1)

    @pytest.mark.e2e
    def test_pwdft_ch4(self):
        """
        Test calculation of CH4
        """
        # create a 10-bohr cubic system with an explicit cutoff
        s = PeriodicSystem(10, ecut=3)
        
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
        np.testing.assert_almost_equal(res['energy'], -18.016869058877973,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-2.87393052, -0.45945895, -0.45945894, -0.45945894,
                -0.22743211]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)
        
    @pytest.mark.e2e
    def test_pwdft_co(self):
        """
        Test calculation of CO
        """
        # create a 10-bohr cubic system with an explicit cutoff
        s = PeriodicSystem(10, ecut=2)
        
        # add methane molecule to system
        atompos = np.array([[5.00000000, 5.00000000, 3.70973478],
                            [5.00000000, 5.00000000, 5.96769844]])
        
        # specify atomic charges
        charges = [6, 8] # C + 4 x H
        
        # add atoms to the system
        for p,c in zip(atompos, charges):
            s.add_atom(p[0], p[1], p[2], c)
        
        # perform SCF calculation
        calculator = PyPWDFT(s, pseudopotential=GTHPseudopotential(s))
        res = calculator.scf(
            tol=1e-4, density_tol=1e-3, nsol=7, verbose=False
        )
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -16.54895159372596,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-1.63632198, -0.61984203, -0.35651651, -0.35649672,
                -0.18413863, 0.02084852, 0.09183999]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)

    @pytest.mark.e2e
    def test_pwdft_co_translation(self):
        """
        Test calculation of CO with translated molecule; the result should be
        - within numerical approximation - the same
        """
        # create a 10-bohr cubic system with an explicit cutoff
        s = PeriodicSystem(10, ecut=2)
        
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
        calculator = PyPWDFT(s, pseudopotential=GTHPseudopotential(s))
        res = calculator.scf(
            tol=1e-4, density_tol=1e-3, nsol=7, verbose=False
        )
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -16.54895159372596,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-1.63632198, -0.61984203, -0.35651651, -0.35649672,
                -0.18413863, 0.02084852, 0.09183999]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)

    @pytest.mark.e2e
    def test_pwdft_co_npfft(self):
        """
        Test calculation of CO using Numpy FFT
        """
        # create a 10-bohr cubic system with an explicit cutoff
        s = PeriodicSystem(10, ecut=2)
        
        # add methane molecule to system
        atompos = np.array([[5.00000000, 5.00000000, 3.70973478],
                            [5.00000000, 5.00000000, 5.96769844]])
        
        # specify atomic charges
        charges = [6, 8] # C + 4 x H
        
        # add atoms to the system
        for p,c in zip(atompos, charges):
            s.add_atom(p[0], p[1], p[2], c)
        
        # perform SCF calculation
        calculator = PyPWDFT(
            s, fft='numpy', pseudopotential=GTHPseudopotential(s)
        )
        res = calculator.scf(
            tol=1e-4, density_tol=1e-3, nsol=7, verbose=False
        )
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -16.54895159372596,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-1.63632198, -0.61984203, -0.35651651, -0.35649672,
                -0.18413863, 0.02084852, 0.09183999]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)
    
    @pytest.mark.e2e
    def test_pwdft_co_scipy(self):
        """
        Test calculation of CO using Numpy FFT
        """
        # create a 10-bohr cubic system with an explicit cutoff
        s = PeriodicSystem(10, ecut=2)
        
        # add methane molecule to system
        atompos = np.array([[5.00000000, 5.00000000, 3.70973478],
                            [5.00000000, 5.00000000, 5.96769844]])
        
        # specify atomic charges
        charges = [6, 8] # C + 4 x H
        
        # add atoms to the system
        for p,c in zip(atompos, charges):
            s.add_atom(p[0], p[1], p[2], c)
        
        # perform SCF calculation
        calculator = PyPWDFT(
            s, fft='scipy', pseudopotential=GTHPseudopotential(s)
        )
        res = calculator.scf(
            tol=1e-4, density_tol=1e-3, nsol=7, verbose=False
        )
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -16.54895159372596,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-1.63632198, -0.61984203, -0.35651651, -0.35649672,
                -0.18413863, 0.02084852, 0.09183999]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)

if __name__ == '__main__':
    unittest.main()
