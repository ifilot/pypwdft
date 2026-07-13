import unittest
import sys
import os
import importlib
import numpy as np
import pytest

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft.gth import GTHPseudopotential
from pypwdft.psystem import PeriodicSystem
from pypwdft.pypwdft import PyPWDFT

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

    def test_initial_arpack_vector_is_flat_for_multiple_bands(self):
        solver_module = importlib.import_module("pypwdft.pypwdft")
        system = PeriodicSystem(10, ecut=1)
        system.add_atom(4.3, 5, 5, 1)
        system.add_atom(5.7, 5, 5, 1)

        class EigensolverReached(Exception):
            pass

        original_eigsh = solver_module.scipy.sparse.linalg.eigsh

        def inspect_eigsh(operator, count, *, which, v0):
            self.assertEqual(v0.shape, (operator.shape[0],))
            self.assertEqual(count, 2)
            raise EigensolverReached

        solver_module.scipy.sparse.linalg.eigsh = inspect_eigsh
        try:
            with self.assertRaises(EigensolverReached):
                PyPWDFT(system, fft="numpy").scf(nsol=2)
        finally:
            solver_module.scipy.sparse.linalg.eigsh = original_eigsh

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
            tol=1e-6, density_tol=1e-6, nsol=7, verbose=False
        )
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -16.54895159372596,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-1.63843684, -0.62127332, -0.35781872, -0.35781871,
                -0.18488378, 0.02060324, 0.09112234]
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
            tol=1e-6, density_tol=1e-6, nsol=7, verbose=False
        )
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -16.54895159372596,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-1.63843684, -0.62127332, -0.35781872, -0.35781871,
                -0.18488378, 0.02060324, 0.09112234]
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
            tol=1e-6, density_tol=1e-6, nsol=7, verbose=False
        )
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -16.54895159372596,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-1.63843684, -0.62127332, -0.35781872, -0.35781871,
                -0.18488378, 0.02060324, 0.09112234]
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
            tol=1e-6, density_tol=1e-6, nsol=7, verbose=False
        )
        
        # test total energy
        np.testing.assert_almost_equal(res['energy'], -16.54895159372596,
                                       decimal=4)
        
        # test eigenvalues
        orbe = [-1.63843684, -0.62127332, -0.35781872, -0.35781871,
                -0.18488378, 0.02060324, 0.09112234]
        np.testing.assert_almost_equal(res['orbe'], orbe, decimal=4)


def test_verbose_scf_prints_calculation_summary(capsys):
    system = PeriodicSystem(10, ecut=1)
    system.add_atom(4.3, 5, 5, 1)
    system.add_atom(5.7, 5, 5, 1)

    with pytest.raises(RuntimeError, match='did not converge'):
        PyPWDFT(system, fft='numpy').scf(maxiter=1, verbose=True)

    output = capsys.readouterr().out
    assert "PyPWDFT self-consistent field calculation" in output
    assert "FFT backend            : NumPy (CPU)" in output
    assert (
        "Cubic unit cell        : 10.000000 x 10.000000 x 10.000000 bohr"
        in output
    )
    assert "Plane waves            :" in output
    assert "Atom positions (bohr)" in output
    assert "Z=1" in output
    assert "Density mixing         : Pulay (fraction=0.5, history=6)" in output
    assert "SCF iterations" in output

if __name__ == '__main__':
    unittest.main()
