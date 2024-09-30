.. _background:
.. index:: Background

Background
==========

Below, a condensed overview is provided on the background of plane wave density
functional theory, mainly scoping these programs with respect to
localized-orbital density functional theory. A detailed treatise on the
mathematical background behind :program:`PyPWDFT` can be found via this link. At
the bottom of this page, useful references are provided to more reading material
and other code sources.

Contemporary electronic structure calculations are predominantly performed using
density functional theory. Two important flavors of density functional theory
(DFT) exist, being localized orbital DFT and plane-wave DFT. In both methods,
the electronic structure problem is solved by means of the introduction of a
basis set. The nature of this basis set however differs.

In localized orbital DFT, the basis functions correspond to quasi atomic
orbitals, such as sets of Slater orbitals or combinations of Gaussian type
orbitals. Furthermore, the atomic system is placed inside an infinite vacuum. In
contrast, in plane wave DFT, it is not the atoms that "spawn" the basis
functions, but the unit cell itself. On top of that, the system is assumed to
be inherently periodic, although "vacuum" or "gas-phase" systems can be
imitated by using relatively large unit cells.

Irrespective of whether one uses localized orbital DFT or plane wave DFT, in the
end the `Roothaan equation <https://en.wikipedia.org/wiki/Roothaan_equations>`_
is being solved. The way this equation is solved however differs. In localized
orbital DFT, the Hamiltonian matrix is fairly small and one can perform its
diagonalization to compute the eigenvalue and -vector pairs. In contrast, in
plane wave DFT, due to the sheer number of plane waves being used, the
Hamiltonian matrix is very large and is in fact not even stored in the
calculation. To find the eigenvector and -value pairs, a partial matrix solver
procedure is used. This matrix solver does not even require the full Hamiltonian
matrix to be known, but only its **effect** on an arbitrary vector
:math:`\vec{c}` in Hilbert space.

Similar to localized-orbital theory, electron-electron interactions can be more
efficiently evaluated by solving the Poisson equation. Since the plane wave
basis set can also be used to describe the electron density, solving for the
Hartree potential effectively utilizes highly-efficient Fast Fourier transforms.
In a similar fashion, also the nuclear attraction energy is solved via
(analytical) Fourier transforms.

Within :program:`PyPWDFT`, all these aspects which make a plane wave DFT code
stand out from a "regular" localized orbital code are shown. Despite the
ambition to show the most salient details, certain aspects are considered out of
scope, though these aspects are important in commercial codes. For example,
:math:`\vec{k}`-point sampling is not shown, nor are pseudopotentials being
discussed. Furthermore, although the code would support non-cubic unit cells,
none of the examples use these. When using the code, please take these
limitations in mind.

Other software
--------------

Besides :program:`PyPWDFT`, there are a number of other educational electronic
structure codes written by this author. Below, a list is provided.

* `PyQInt <https://pyqint.imc-tue.nl/>`_ : Hartree-Fock code with many additional features.
* `PyDFT <https://pydft.imc-tue.nl/>`_ : Localized orbital density functional theory code.
* `HFCXX <https://github.com/ifilot/hfcxx>`_ : Hartree-Fock code written in C++.
* `DFTCXX <https://github.com/ifilot/dftcxx>`_ : Localized-orbital density functional theory code written in C++.

Further reading
---------------

* A detailed background on electronic structure calculations, including a
  derivation of the Hartree-Fock and localized-orbital density functional theory
  calculations can be found in `this open-access textbook
  <https://ifilot.pages.tue.nl/elements-of-electronic-structure-theory/>`_.