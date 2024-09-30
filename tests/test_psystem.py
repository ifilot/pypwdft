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
        # create cubic periodic system with lattice size of 10 A and
        # 32 grid points per cartesian direction
        s = PeriodicSystem(10, 32)
        
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
        # create cubic periodic system with lattice size of 10 A and
        # 32 grid points per cartesian direction
        s = PeriodicSystem(10, 32)
        
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

if __name__ == '__main__':
    unittest.main()
