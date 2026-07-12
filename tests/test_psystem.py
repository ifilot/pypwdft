import unittest
import sys
import os
import numpy as np

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import PeriodicSystem

class TestPeriodicUnitCell(unittest.TestCase):
        
    def test_unitcell_properties(self):
        """
        Test whether the PeriodicUnitCell class has correctly stored
        all relevant properties that describe the unit cell
        """
        # create a 10-bohr cubic system with a 5-Hartree cutoff
        s = PeriodicSystem(10, ecut=5)
        
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
        
        # test size of unit cell
        np.testing.assert_equal(s.get_omega(), 1000)
        
        # test that all atoms are entered correctly in the system
        np.testing.assert_almost_equal(s.get_atom_positions(), atompos)
        
        # check the charges
        np.testing.assert_almost_equal(s.get_atom_charges(), charges)
        
        # check the size of the k2 array
        npts = s.get_npts()
        np.testing.assert_almost_equal(s.get_pw_k2().shape, [npts, npts, npts])

    def test_translation(self):
        s = PeriodicSystem(10, ecut=5)
        
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

        s.translate((5,5,5))

        # because of modulus, translation by (5,5,5) should yield (0,0,0)
        np.testing.assert_almost_equal(s.get_atom_positions()[0],
                                       np.array([0.00000000, 0.00000000, 0.00000000]))
        
        # test for fourth H atom
        np.testing.assert_almost_equal(s.get_atom_positions()[-1],
                                       np.array([1.19575624, 8.80424376, 8.80424376]))

    def test_spherical_plane_wave_cutoff(self):
        s = PeriodicSystem(10, ecut=5)

        mask = s.get_pw_mask()
        np.testing.assert_array_less(
            0.5 * s.get_pw_k2()[mask], 5 + np.finfo(float).eps
        )
        self.assertEqual(s.get_n_plane_waves(), np.count_nonzero(mask))
        self.assertLess(s.get_n_plane_waves(), s.get_npts()**3)
        self.assertEqual(s.get_density_ecut(), 20)
        self.assertEqual(s.get_npts(), s.get_density_npts())
        self.assertLess(s.get_wavefunction_npts(), s.get_density_npts())

        finer = PeriodicSystem(10, ecut=5, density_ecut=30)
        self.assertGreater(finer.get_density_npts(), s.get_density_npts())
        self.assertEqual(finer.get_n_plane_waves(), s.get_n_plane_waves())

        with self.assertRaises(ValueError):
            PeriodicSystem(10, ecut=5, density_ecut=19)
        with self.assertRaises(ValueError):
            PeriodicSystem(10, ecut=0)
        with self.assertRaises(TypeError):
            PeriodicSystem(10)
        with self.assertRaises(TypeError):
            PeriodicSystem(10, 16)

if __name__ == '__main__':
    unittest.main()
