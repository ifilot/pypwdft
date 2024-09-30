import unittest
import sys
import os
import numpy as np

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import PeriodicSystem, PyPWDFT, SystemBuilder

class TestSystemBuilder(unittest.TestCase):
    """
    Test evaluation of terms in reciprocal space
    """

    def test_system_builder(self):
        """
        Test equivalence assessment real and reciprocal space Hartree energy
        """
        # Grab system from SystemBuilder
        s = SystemBuilder().from_name('CH4', sz = 10.0, npts = 32)

        atompos = np.array([[5.00000000, 5.00000000, 5.00000000],
                            [6.19575624, 6.19575624, 6.19575624],
                            [3.80424376, 3.80424376, 6.19575624],
                            [3.80424376, 6.19575624, 3.80424376],
                            [6.19575624, 3.80424376, 3.80424376]])
        
        atomcharges = [6,1,1,1,1]

        np.testing.assert_almost_equal(s.get_atom_positions(), atompos, decimal=4)
        np.testing.assert_almost_equal(s.get_atom_charges(), atomcharges, decimal=4)

if __name__ == '__main__':
    unittest.main()
