.. _usage:
.. index:: Usage

Usage
=====

Creating atomic systems
-----------------------

To perform an electronic structure calculation, one first has to define an
atomic system. This can be done in two ways, either manually or by means of the
:class:`pypwdft.SystemBuilder` class.

Automatic
#########

The most straightforward way is to use the :class:`pypwdft.SystemBuilder` class.
To build a system, one can run

.. code:: python

    from pypwdft import SystemBuilder, PeriodicSystem

    # create cubic periodic system with lattice size of 10 Bohr units
    npts = 16   # number of grid points
    sz = 10
    
    # construct CH4 molecule system via SystemBuilder
    s = SystemBuilder().from_name('CH4', sz=sz, npts=npts)

To view the unit cell matrix and the atomic coordinates, one can simply invoke

.. code:: python

    print(s)

which in the above scenario yields the following output

.. code::

    [[10.  0.  0.]
     [ 0. 10.  0.]
     [ 0.  0. 10.]]
    6  (5.0000  5.0000  5.0000)
    1  (6.1958  6.1958  6.1958)
    1  (3.8042  3.8042  6.1958)
    1  (3.8042  6.1958  3.8042)
    1  (6.1958  3.8042  3.8042)

The following molecules are available via the :class:`pypwdft.SystemBuilder`
class:

* benzene
* bf3
* ch4
* co
* ethylene
* h2
* h2o
* he
* lih
* nh3

.. note::

    Unless otherwise specified, atomic units are used throughout the program.
    This means that all distances are provided in Bohr units.

Manual
######

Alternatively, one can also build a system by hand. First, define the unit cell
and the number of sampling points per Cartesian direction.

.. code:: python

    from pypwdft import PeriodicSystem

    npts = 32   # number of grid points
    sz = 10     # edge size of cubic unit cell
    s = PeriodicSystem(sz, npts)

Next, one can add atoms to the PeriodicSystem by means of the 
:class:`pypwdft.PeriodicSystem.add_atom` method

.. code:: python

    # atomic positions
    atompos = np.array([[5.00000000, 5.00000000, 5.00000000],
                        [6.19575624, 6.19575624, 6.19575624],
                        [3.80424376, 3.80424376, 6.19575624],
                        [3.80424376, 6.19575624, 3.80424376],
                        [6.19575624, 3.80424376, 3.80424376]])
    
    # atomic charges
    charges = [6, 1, 1, 1, 1] # C + 4 x H

    # add atoms to the system
    for p,c in zip(atompos, charges):
        s.add_atom(p[0], p[1], p[2], c)

Performing electronic structure calculation
-------------------------------------------

Electronic structure calculations are handled by the :class:`pypwdft.PyPWDFT`
class. For each separate electronic structure calculation, one creates a fresh
instance of this class. Upon instancing, a :class:`pypwdft.PeriodicSystem`
instance is provided. Furthermore, the user can select which FFT algorithm is
being used. Three options are available:

* `NumPy FFT <https://numpy.org/doc/stable/reference/routines.fft.html>`_ : :code:`numpy`
* `Scipy FFT <https://docs.scipy.org/doc/scipy/tutorial/fft.html>`_ : :code:`scipy`
* `pyFFTW FFT <https://pypi.org/project/pyFFTW/>`_ : :code:`pyfftw`

By default, :code:`pyfftw` is used as this algorithm provides the best
performance.

After building the :class:`pypwdft.PyPWDFT` instance, the electronic structure
calculation can be invoked by running :class:`pypwdft.PyPWDFT.scf`. This
function allows for several parameters tuning the execution. The default
parameters are however typically suitable. By default, no output is written to
the console, but this can be changed by setting :code:`verbose=True`. Below,
an example is provided how to set-up a electronic structure calculation.

.. code:: python

    from pypwdft import PyPWDFT, PeriodicSystem, SystemBuilder

    # create cubic periodic system with lattice size of 10 Bohr units
    npts = 32   # number of grid points
    sz = 10
    
    # construct CH4 molecule system via SystemBuilder
    s = SystemBuilder().from_name('CH4', sz=sz, npts=npts)
        
    # construct calculator object
    calculator = PyPWDFT(s)
    
    # perform self-consistent field procedure and store results in res object
    res = calculator.scf(tol=1e-5, verbose=True)

Upon execution of this code, the following output is written to the console

.. code::

    001 | Etot =  25.72188838 Ht | eps = 2.5722e+01 | dt = 0.5023 s
    002 | Etot =  -3.73540233 Ht | eps = 2.9457e+01 | dt = 0.5098 s
    003 | Etot = -21.73526459 Ht | eps = 1.8000e+01 | dt = 0.5807 s
    004 | Etot = -29.42671403 Ht | eps = 7.6914e+00 | dt = 0.7262 s
    005 | Etot = -32.91894875 Ht | eps = 3.4922e+00 | dt = 0.8088 s
    006 | Etot = -34.79324803 Ht | eps = 1.8743e+00 | dt = 0.8830 s
    007 | Etot = -35.89028706 Ht | eps = 1.0970e+00 | dt = 0.9239 s
    008 | Etot = -36.56316890 Ht | eps = 6.7288e-01 | dt = 0.9665 s
    009 | Etot = -36.98653043 Ht | eps = 4.2336e-01 | dt = 1.0637 s
    010 | Etot = -37.25642182 Ht | eps = 2.6989e-01 | dt = 0.9757 s
    011 | Etot = -37.42963747 Ht | eps = 1.7322e-01 | dt = 0.8910 s
    012 | Etot = -37.54120394 Ht | eps = 1.1157e-01 | dt = 0.9177 s
    013 | Etot = -37.61320317 Ht | eps = 7.1999e-02 | dt = 0.8708 s
    014 | Etot = -37.65971654 Ht | eps = 4.6513e-02 | dt = 0.9530 s
    015 | Etot = -37.68978030 Ht | eps = 3.0064e-02 | dt = 0.9666 s
    016 | Etot = -37.70921444 Ht | eps = 1.9434e-02 | dt = 1.0423 s
    017 | Etot = -37.72177570 Ht | eps = 1.2561e-02 | dt = 1.0124 s
    018 | Etot = -37.72989204 Ht | eps = 8.1163e-03 | dt = 0.9978 s
    019 | Etot = -37.73513379 Ht | eps = 5.2417e-03 | dt = 1.1720 s
    020 | Etot = -37.73851690 Ht | eps = 3.3831e-03 | dt = 1.2111 s
    021 | Etot = -37.74069871 Ht | eps = 2.1818e-03 | dt = 1.5268 s
    022 | Etot = -37.74210447 Ht | eps = 1.4058e-03 | dt = 1.0254 s
    023 | Etot = -37.74300923 Ht | eps = 9.0476e-04 | dt = 1.0217 s
    024 | Etot = -37.74359080 Ht | eps = 5.8157e-04 | dt = 1.0133 s
    025 | Etot = -37.74396409 Ht | eps = 3.7329e-04 | dt = 1.3798 s
    026 | Etot = -37.74420328 Ht | eps = 2.3920e-04 | dt = 1.2190 s
    027 | Etot = -37.74435626 Ht | eps = 1.5298e-04 | dt = 1.2910 s
    028 | Etot = -37.74445388 Ht | eps = 9.7624e-05 | dt = 1.0146 s
    029 | Etot = -37.74451603 Ht | eps = 6.2144e-05 | dt = 1.0786 s
    030 | Etot = -37.74455548 Ht | eps = 3.9447e-05 | dt = 0.8716 s
    031 | Etot = -37.74458043 Ht | eps = 2.4958e-05 | dt = 0.9166 s
    032 | Etot = -37.74459617 Ht | eps = 1.5732e-05 | dt = 1.0046 s
    033 | Etot = -37.74460604 Ht | eps = 9.8741e-06 | dt = 1.0156 s

Furthermore, the result of the calculation are stored in a so-called results
dictionary. This dictionary contains the following entries:

* :code:`Energy` : List of total electronic energy per iteration
* :code:`Etot` : Total electronic energy
* :code:`Ekin` : Kinetic energy
* :code:`Enuc` : Nuclear attraction energy
* :code:`Erep` : Electron-electron repulsion energy
* :code:`Exc` : Exchange-correlation energy
* :code:`edens` : Electron density scalar field
* :code:`k2` : Plane wave vector magnitudes
* :code:`dV` : Real-space integration constant
* :code:`Eewald` : Nuclear-nuclear repulsion (Ewald sum)
* :code:`orbc_fft` : Reciprocal-space representation of the molecular orbitals
* :code:`orbe` : Molecular orbital energies
* :code:`orbc_rs` : Real-space representation of the molecular orbitals
* :code:`ttime` : Total computation time

.. note::

    All energies are provided in Hartrees.

Visualizing molecular orbitals
------------------------------

Molecular orbitals can be visualized using `matplotlib
<https://matplotlib.org/>`_. The results of the :class:`pypwdft.PyPWDFT.scf`
routine is a dictionary of which one of its elements is :code:`orbc_rs`. This
element corresponds to a four-dimensional array where the first index loops over
the molecular orbitals with increasing orbital energy. 

Two-dimensional projections
###########################

By setting the two indices, one can essentially extract a specific z-layer from
the scalar field of any of the molecular orbitals and in turn visualize these
using the :code:`imshow` function of matplotlib. This is demonstrated in the
script below.

.. note::

    The wave functions generated by a plane wave density functional theory
    calculation are complex values, i.e.

    .. math::

        \psi(\vec{r}) \in \mathbb{C}^{3}

    and as such, it is recommended to visualize both the real and complex parts of
    the wave function as well as its electron density as given by

    .. math::

        \rho = |\psi(\vec{r})|^{2} = \psi(\vec{r}) \cdot \psi^{*}(\vec{r}).

.. code::

    # import the required libraries for the test
    from pypwdft import PyPWDFT, PeriodicSystem, SystemBuilder
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    def main():
        # create cubic periodic system with lattice size of 10 Bohr units
        npts = 32   # number of grid points
        sz = 10
        
        # construct CH4 molecule system via SystemBuilder
        s = SystemBuilder().from_name('CH4', sz=sz, npts=npts)
            
        # construct calculator object
        calculator = PyPWDFT(s)
        
        # perform self-consistent field procedure and store results in res object
        res = calculator.scf(tol=1e-1, verbose=True)
        
        # visualize the occupied molecular orbitals
        fig, im = plt.subplots(3,5, dpi=300, figsize=(16,8))
        extent=[-sz/2,sz/2,-sz/2,sz/2]
        orbe = res['orbe']
        fig.suptitle('Occupied molecular orbitals of CH$_{4}$')
        m = np.empty((3,5), dtype=object) # create placeholder for maps
        
        for i in range(0,5):
            # visualize the real part of the wave function
            field = np.real(res['orbc_rs'][i][npts//2, :, :])
            maxval = max(np.max(np.abs(field)), 0.1)
            m[0][i] = im[0,i].imshow(field, origin='lower',
                    interpolation='bicubic', extent=extent, cmap='PiYG',
                    vmin=-maxval, vmax=maxval)
            im[0,i].set_title(r'$\mathbb{R}\;[\psi_{%i}]$' % (i+1))
            
            # visualize the imaginary part of the wave function
            field = np.imag(res['orbc_rs'][i][npts//2, :, :])
            maxval = max(np.max(np.abs(field)), 0.001)
            m[1][i] = im[1,i].imshow(field, origin='lower',
                    interpolation='bicubic', extent=extent, cmap='PiYG',
                    vmin=-maxval, vmax=maxval)
            im[1,i].set_title(r'$\mathbb{I}\;[\psi_{%i}]$' % (i+1))
            
            # visualize the electron density
            m[2][i] = im[2,i].imshow(np.real(res['orbc_rs'][i][npts//2, :, :].conj() * 
                                res['orbc_rs'][i][npts//2, :, :]), 
                        origin='lower', interpolation='bicubic', extent=extent)
            im[2,i].set_title(r'$\rho_{%i}$' % (i+1))
            
        for j in range(0,3):
            for i in range(0,5):
                im[j,i].set_xlabel('$x$-position [a.u.]')
                im[j,i].set_ylabel('$y$-position [a.u.]')
                
                divider = make_axes_locatable(im[j,i])
                cax = divider.append_axes('right', size='5%', pad=0.05)
                fig.colorbar(m[j][i], cax=cax, orientation='vertical')
            
        plt.tight_layout()
        
    if __name__ == '__main__':
        main()

The result of running this script is shown below.

.. image:: _static/img/orbs_ch4.png

Isosurfaces
###########

Alternative to two-dimensional projections, one can also create isosurfaces. For
this, we will use the external module `PyTessel
<https://pytessel.imc-tue.nl/>`_. For isosurfaces, a relatively large number of
sampling points for the scalar fields are required. However, this comes at the
expense of computational time and for this reason, we will perform the
electronic structure calculation initially using only 32 sampling points per
Cartesian direction and follow up using quintic interpolation to "upsample" the
scalar fields. An example of this process is shown in the image below.

.. code::

    import numpy as np
    from pytessel import PyTessel
    from scipy.interpolate import RegularGridInterpolator
    from pypwdft import PyPWDFT, SystemBuilder, PeriodicSystem

    def main():
        # create cubic periodic system with lattice size of 10 Bohr units
        npts = 32       # number of grid points
        sz = 10
        
        # construct CO molecule system via SystemBuilder
        s = SystemBuilder().from_name('CO', sz=sz, npts=npts)
            
        # construct calculator object
        calculator = PyPWDFT(s)
        
        # perform self-consistent field procedure and store results in res object
        res = calculator.scf(tol=1e-5, nsol=9, verbose=True)
        
        # generate PyTessel object
        pytessel = PyTessel()
        
        for i in range(2,9):
            print('Building isosurfaces: %02i' % (i+1))
            scalarfield = interpolate_grid(res['orbc_rs'][i], sz, npts, 4)
            unitcell = np.identity(3) * sz
            
            # build positive real isosurface
            vertices, normals, indices = pytessel.marching_cubes(scalarfield.real.flatten(), scalarfield.shape, unitcell.flatten(), 0.1)
            pytessel.write_ply('MO_PR_%02i.ply' % (i+1), vertices, normals, indices)
            
            # build negative real isosurface
            vertices, normals, indices = pytessel.marching_cubes(scalarfield.real.flatten(), scalarfield.shape, unitcell.flatten(), -0.1)
            pytessel.write_ply('MO_NR_%02i.ply' % (i+1), vertices, normals, indices)
            
            # build positive imaginary isosurface
            vertices, normals, indices = pytessel.marching_cubes(scalarfield.imag.flatten(), scalarfield.shape, unitcell.flatten(), 0.1)
            pytessel.write_ply('MO_PI_%02i.ply' % (i+1), vertices, normals, indices)
            
            # build negative imaginary isosurface
            vertices, normals, indices = pytessel.marching_cubes(scalarfield.imag.flatten(), scalarfield.shape, unitcell.flatten(), -0.1)
            pytessel.write_ply('MO_NI_%02i.ply' % (i+1), vertices, normals, indices)

    def interpolate_grid(scalarfield, sz, npts, amp=2):
        x = np.linspace(0, sz, npts)
        interp = RegularGridInterpolator((x,x,x), scalarfield, method='quintic')
        s = PeriodicSystem(sz, npts * amp)
        points = s.get_r()
        
        return interp(points)

    if __name__ == '__main__':
        main()

In the image below, the isosurfaces corresponding to the real part of the scalar
field for the 3σ, 4σ and 1π orbitals are visualized.

.. image:: _static/img/orbs_co_isosurfaces.png