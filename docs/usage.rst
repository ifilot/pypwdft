.. _usage:
.. index:: Usage

Usage
=====

Quick start
-----------

A calculation consists of a :class:`pypwdft.Structure`, a
:class:`pypwdft.PWDFT` configuration, and the resulting
:class:`pypwdft.DFTResult`:

.. code:: python

    from pypwdft import PWDFT, Structure

    structure = Structure.from_name("H2", cell=10)
    calculation = PWDFT(
        structure,
        cutoff=40,
        xc="pbe",
        pseudopotential="gth",
    )
    result = calculation.run(convergence=1e-6, verbosity=1)

    print(result.energy.total)
    result.plot_orbitals(plane="xz", save="h2-orbitals.pdf")

PyPWDFT uses atomic units: distances are in bohr and energies are in Hartree
unless another unit is explicitly selected.

Structures
----------

Bundled molecules can be loaded by name. They are centred automatically in a
cubic cell:

.. code:: python

    methane = Structure.from_name("CH4", cell=12)

Available bundled structures include ``benzene``, ``bf3``, ``bh3``, ``ch4``,
``co``, ``co2``, ``ethylene``, ``h2``, ``h2o``, ``he``, ``lih``, and ``nh3``.

XYZ coordinates are interpreted as angstrom while ``cell`` remains expressed
in bohr:

.. code:: python

    water = Structure.from_xyz("water.xyz", cell=12)

A structure can also be constructed directly. Setting ``units="angstrom"``
applies to both the coordinates and cell:

.. code:: python

    hydrogen = Structure(
        symbols=["H", "H"],
        positions=[[0, 0, -0.37], [0, 0, 0.37]],
        cell=6,
        units="angstrom",
    )

Calculation settings
--------------------

The wavefunction cutoff and optional density cutoff belong to the calculation,
not the geometry:

.. code:: python

    calculation = PWDFT(
        methane,
        cutoff=40,
        density_cutoff=160,
        xc="pbe",
        pseudopotential="gth",
        device="cpu",
    )

``xc`` accepts ``"lda"`` and ``"pbe"``. The ``"gth"`` pseudopotential
selection automatically uses GTH-PADE for LDA and GTH-PBE for PBE. An explicit
configuration supports custom parameter directories and valence choices:

.. code:: python

    from pypwdft import GTH

    calculation = PWDFT(
        structure,
        cutoff=40,
        xc="pbe",
        pseudopotential=GTH(charges={"Ga": 13}),
    )

Use ``pseudopotential="all-electron"`` to disable the frozen-core model.

SCF settings
------------

Common settings can be passed directly to :meth:`pypwdft.PWDFT.run`:

.. code:: python

    result = calculation.run(
        convergence=1e-6,
        density_convergence=1e-6,
        max_iterations=150,
        bands=8,
        verbosity=1,
    )

Reusable settings can be collected in :class:`pypwdft.SCFSettings`:

.. code:: python

    from pypwdft import SCFSettings

    settings = SCFSettings(
        convergence=1e-6,
        density_convergence=1e-6,
        max_iterations=150,
        bands=8,
        verbosity=1,
    )
    result = calculation.run(settings)

Results
-------

Results expose named, typed groups instead of abbreviated dictionary keys:

.. code:: python

    result.energy.total
    result.energy.kinetic
    result.energy.electron_ion
    result.energy.nonlocal_
    result.energy.hartree
    result.energy.xc
    result.energy.ion_ion

    result.density
    result.orbitals.energies
    result.orbitals.real_space
    result.orbitals.reciprocal_space

    result.basis.cutoff
    result.basis.density_cutoff
    result.basis.plane_waves
    result.basis.wavefunction_grid
    result.basis.density_grid

    result.scf.converged
    result.scf.iterations
    result.scf.energy_residual
    result.scf.density_residual
    result.scf.elapsed_time
    result.scf.backend

Plotting
--------

Occupied orbitals are plotted by default:

.. code:: python

    result.plot_orbitals(
        plane="xz",
        columns=3,
        save="orbitals.pdf",
        show_imaginary_norm=False,
    )

Residual imaginary norms are shown in orbital titles by default. Setting
``show_imaginary_norm=False`` hides that diagnostic without changing the
orbital or the plotted field.

For a collection containing different nodal planes, ``plane="auto"`` selects
the central Cartesian plane with the largest norm separately for each orbital.
Use ``ngrid`` and ``tick_rotation`` to keep tick labels readable in dense
subplot grids.

Specific orbitals and the converged density can also be selected:

.. code:: python

    result.plot_orbitals(indices=[0, 2, 4], plane="xy")
    result.plot_density(plane="xz", save="density.pdf")

GPU calculations
----------------

Install the CuPy package compatible with the local CUDA runtime, then request
the CUDA device explicitly:

.. code:: python

    calculation = PWDFT(
        structure,
        cutoff=40,
        xc="pbe",
        pseudopotential="gth",
        device="cuda",
    )

Returned arrays are NumPy arrays, independent of the calculation device.
