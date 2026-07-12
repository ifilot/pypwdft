import unittest
import sys
import os
import numpy as np

# add a reference to load the PyPWDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# import the required libraries for the test
from pypwdft import SystemBuilder

class TestPeriodicUnitCell(unittest.TestCase):

    def test_ewald_CH4(self):
        """
        Test calculation of Ewald sum
        """
        # The Ewald energy is independent of the automatically selected grid.
        s = SystemBuilder().from_name('CH4', sz=10, ecut=1)
        
        # calculate Ewald sum
        Eewald = s.calculate_ewald_sum()
        
        # test that the Ewald sum is correctly calculated
        np.testing.assert_almost_equal(Eewald, -0.480804464127026)

    def test_ewald_CO(self):
        """
        Test calculation of Ewald sum
        """
        # The Ewald energy is independent of the automatically selected grid.
        s = SystemBuilder().from_name('CO', sz=10, ecut=1)
        
        # calculate Ewald sum
        Eewald = s.calculate_ewald_sum()
        
        # test that the Ewald sum is correctly calculated
        np.testing.assert_almost_equal(Eewald, -5.995622350717917)
        
if __name__ == '__main__':
    unittest.main()
