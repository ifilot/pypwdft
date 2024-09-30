# -*- coding: utf-8 -*-

# 
# This file is part of the PyPWDFT distribution 
# Copyright (c) 2024 Ivo Filot
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import os
from .psystem import PeriodicSystem
from mendeleev import element

class SystemBuilder:
    """
    Class that builds molecules from templates or from point group descriptions
    """
    def __init__(self):
        """Create the object
        """
        pass

    def from_name(self, molname:str, sz:float=10, npts:int=32) -> PeriodicSystem:
        """Create a PeriodicSystem instance by specifying the name of the molecule

        Args:
            molname (str): name of the molecule
            sz (float, optional): edge size of the cubic unit cell in a.u. Defaults to 10 a.u.
            npts (int, optional): number of sampling points per Cartesian direction. Defaults to 32.

        Returns:
            PeriodicSystem: PeriodicSystem instance encapsulating the specified molecule
        """
        
        fname = os.path.join(os.path.dirname(__file__), 'molecules', molname.lower() + '.xyz')
        return self.from_file(fname, sz=sz, npts=npts)
       
    def from_file(self, path:str, sz:float=10, npts:int=32) -> PeriodicSystem:
        """Construct a PeriodicSystem instance from a .xyz file

        Args:
            path (str): path to .xyz file
            sz (float, optional): edge size of the cubic unit cell in a.u. Defaults to 10 a.u.
            npts (int, optional): number of sampling points per Cartesian direction. Defaults to 32.

        Returns:
            PeriodicSystem: PeriodicSystem instance encapsulating the specified molecule
        """
        with open(path, 'r') as f:
            lines = f.readlines()
            nratoms = int(lines[0].strip())
            psys = PeriodicSystem(sz=sz, npts=npts)
            hsz = sz * 0.5

            for line in lines[2:2+nratoms]:
                pieces = line.split()

                el = element(pieces[0])

                psys.add_atom(hsz + float(pieces[1]) * 1.8897259886, 
                              hsz + float(pieces[2]) * 1.8897259886, 
                              hsz + float(pieces[3]) * 1.8897259886, 
                              charge=el.atomic_number, unit='bohr')

            return psys
