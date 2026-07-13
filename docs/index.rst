PyPWDFT: pure-python plane-wave density functional theory
=========================================================

.. image:: https://anaconda.org/ifilot/pypwdft/badges/version.svg
   :target: https://anaconda.org/ifilot/pypwdft
.. image:: https://img.shields.io/pypi/v/pypwdft?color=green
   :target: https://pypi.org/project/pypwdft/
.. image:: https://github.com/ifilot/pypwdft/actions/workflows/build_pypi.yml/badge.svg
   :target: https://github.com/ifilot/pypwdft/actions/workflows/build_pypi.yml
.. image:: https://github.com/ifilot/pypwdft/actions/workflows/build_conda.yml/badge.svg
   :target: https://github.com/ifilot/pypwdft/actions/workflows/build_conda.yml
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

:program:`PyPWDFT` is a pure-Python package for performing plane-wave DFT
calculations. :program:`PyPWDFT` is designed with the aim of teaching students
the inner workings of a PW-DFT code. The source code is kept relatively small
and contains extensive commenting. Using :program:`PyPWDFT` is fairly simple,
though its application is due to its scope relatively limited.

Features and scope
------------------

* Unit cells are limited to cubes.
* GTH frozen-core pseudopotentials with local and non-local terms.
* Sampling is currently limited to the :math:`\Gamma`-point.
* Lowest Kohn-Sham states are found using the Arnoldi iterative procedure.
* Slater exchange functional
* Vosko-Wilk-Nusair correlation functional (VWN5)
* Perdew-Burke-Ernzerhof (PBE) exchange-correlation functional
* Dualism: the same basis set is used to describe both the molecular orbitals
  as well as the electron density.
* :class:`pypwdft.Structure` loads bundled molecules and XYZ structures.
* The self-consistent field procedure allows for verbose (detailed output)
* Variables and matrices relevant to the computation are accessible to the user
  such that they can follow the procedure.

Example
-------

The code below will perform a plane-wave density functional theory calculation
for the methane molecule inside a 10x10x10 a.u. unit cell.

.. code:: python

   from pypwdft import PWDFT, Structure

   def main():
      structure = Structure.from_name("CH4", cell=10)
      calculation = PWDFT(
         structure,
         cutoff=5,
         xc="pbe",
         pseudopotential="gth",
      )
      result = calculation.run(convergence=1e-4, verbosity=1)
      result.plot_orbitals(plane="xz", save="ch4-orbitals.pdf")

   if __name__ == '__main__':
      main()

The set of molecular orbitals obtained from a methane calculation is shown
below.

.. image:: _static/img/orbs_ch4.png

:program:`PyPWDFT` has been developed at the Eindhoven University of Technology,
Netherlands. :program:`PyPWDFT` and its development are hosted on `Github
<https://www.github.com/ifilot/pypwdft>`_.  Bugs and feature
requests are ideally submitted via the `github issue tracker
<https://www.github.com/ifilot/pypwdft/issues>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   background
   usage
   api
   community_guidelines

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
