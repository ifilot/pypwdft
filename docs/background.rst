.. _background:
.. index:: Background

Background
==========

PyPWDFT is an educational implementation of spin-unpolarized, plane-wave
Kohn--Sham density functional theory. This chapter outlines the equations and
numerical choices used by the program. It is intended to connect the source
code to the underlying theory; it is not a general review of every technique
used by production solid-state codes.

Periodic plane-wave representation
----------------------------------

A plane-wave calculation describes a periodically repeated unit cell. PyPWDFT
currently uses a cubic cell of volume :math:`\Omega=L^3` and samples only the
:math:`\Gamma` point. Molecules are represented by placing them in a cell with
enough vacuum that interactions between periodic images become small. Cell-size
convergence should therefore be checked just as cutoff convergence is checked.

The normalized plane waves are

.. math::

   \langle \mathbf r | \mathbf G \rangle =
   \frac{1}{\sqrt{\Omega}}e^{i\mathbf G\cdot\mathbf r},

where :math:`\mathbf G` is a reciprocal lattice vector. An orbital is expanded
as

.. math::

   \psi_n(\mathbf r) = \frac{1}{\sqrt{\Omega}}
   \sum_{\mathbf G} c_{n\mathbf G}e^{i\mathbf G\cdot\mathbf r}.

Unlike a localized atomic-orbital basis, the number of basis functions is
determined by the cell and a kinetic-energy cutoff. PyPWDFT retains the
spherical set

.. math::

   \frac{1}{2}|\mathbf G|^2 \leq E_\mathrm{cut}.

Increasing :math:`E_\mathrm{cut}` systematically enlarges the basis. It also
increases the FFT grids and computational cost, so energies and other reported
properties should be tested for cutoff convergence.

Wavefunction and density grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Products of wavefunctions contain reciprocal components up to twice the
largest wavefunction vector. Their kinetic-energy scale is consequently four
times larger. The density cutoff therefore defaults to

.. math::

   E_\mathrm{cut}^{\rho} = 4E_\mathrm{cut}.

PyPWDFT derives FFT-friendly wavefunction and density grid sizes from these two
cutoffs. The current solver tracks both sizes but uses the density grid as its
common working FFT grid. This avoids aliasing orbital products while retaining
the spherical :math:`E_\mathrm{cut}` mask as the actual variational orbital
basis.

Kohn--Sham problem
------------------

For a closed-shell system, PyPWDFT solves

.. math::

   \hat H_\mathrm{KS}\psi_n = \varepsilon_n\psi_n,

with

.. math::

   \hat H_\mathrm{KS} =
   -\frac{1}{2}\nabla^2
   + v_\mathrm{ion}^{\mathrm{local}}(\mathbf r)
   + v_\mathrm{H}[\rho](\mathbf r)
   + v_\mathrm{xc}[\rho](\mathbf r)
   + \hat V_\mathrm{ion}^{\mathrm{nonlocal}}.

The density is constructed from the occupied orbitals as

.. math::

   \rho(\mathbf r) = 2\sum_{n=1}^{N_\mathrm{occ}}
   |\psi_n(\mathbf r)|^2.

The factor two accounts for paired spins. The present implementation therefore
requires an even number of electrons and does not treat spin polarization,
fractional occupations, or finite-temperature smearing.

The Hamiltonian is not assembled as a dense matrix. Instead, the code defines
its action on a vector of active plane-wave coefficients. The kinetic term is
diagonal in reciprocal space, while local potentials are applied by
transforming to real space, multiplying by the potential, and transforming
back. The selected CPU or GPU backend's iterative eigensolver obtains only the
lowest requested states. Additional unoccupied states can be requested through
the ``bands`` SCF setting, but only occupied states contribute to the density.

Hartree potential
-----------------

The classical electrostatic potential of the electrons follows from Poisson's
equation,

.. math::

   \nabla^2 v_\mathrm{H}(\mathbf r) = -4\pi\rho(\mathbf r).

In reciprocal space this becomes, for :math:`\mathbf G\neq 0`,

.. math::

   v_\mathrm{H}(\mathbf G) =
   \frac{4\pi}{|\mathbf G|^2}\rho(\mathbf G).

The reciprocal zero component is set to zero, corresponding to the usual
periodic electrostatic reference. FFTs make the transformations between
:math:`\rho(\mathbf r)` and :math:`\rho(\mathbf G)` efficient.

Ions and GTH pseudopotentials
-----------------------------

The recommended ionic model is the Goedecker--Teter--Hutter (GTH)
pseudopotential. Core electrons are replaced by an effective potential and
only the selected valence electrons enter the Kohn--Sham calculation. A GTH
pseudopotential contains

* a local ionic potential, evaluated on the reciprocal grid;
* separable non-local projectors, applied directly to the active plane-wave
  coefficients; and
* valence ionic charges used for the periodic ion--ion energy.

Bundled GTH-PADE parameters are paired with LDA, while bundled GTH-PBE
parameters are paired with PBE. The high-level API enforces this matching unless
an explicit external parameter directory is supplied.

An all-electron periodic Coulomb model remains available for educational and
testing purposes. In that mode the electron--nuclear potential is evaluated in
reciprocal space and the nuclear repulsion is calculated using an Ewald sum.
Core states vary rapidly near a nucleus, however, and require much higher
plane-wave cutoffs; GTH calculations are therefore substantially more practical
for elements heavier than hydrogen and helium.

Exchange and correlation
------------------------

PyPWDFT provides two spin-unpolarized exchange-correlation choices:

``lda``
   Slater exchange combined with the Vosko--Wilk--Nusair VWN5
   parameterization of Ceperley--Alder correlation data.

``pbe``
   The Perdew--Burke--Ernzerhof generalized-gradient approximation. PBE depends
   on both the density and its gradient; derivatives are evaluated using the
   periodic reciprocal representation.

The exchange-correlation energy and potential are evaluated consistently for
the selected functional. PBE calculations should use GTH-PBE parameters, and
LDA calculations should use GTH-PADE parameters.

Self-consistent field procedure
-------------------------------

Because the Hartree and exchange-correlation potentials depend on the density,
the eigenvalue problem must be solved self-consistently. PyPWDFT uses the
following cycle:

#. Start from a homogeneous density with the correct electron count.
#. Construct the Hartree and exchange-correlation potentials.
#. Build the matrix-free Kohn--Sham Hamiltonian.
#. Solve for the requested lowest Kohn--Sham states.
#. Construct a new density from the occupied states.
#. Form the density residual and use Pulay mixing to minimize a short history
   of recent residuals. The first update and ill-conditioned extrapolations
   fall back to damped linear mixing.
#. Repeat until both the total-energy change and RMS density residual satisfy
   their convergence thresholds.

The default Pulay mixer retains six densities and applies a mixing fraction of
0.5. This typically suppresses oscillations and reaches self-consistency in
fewer iterations than a fixed linear update. ``mixing="linear"`` remains
available for comparison and troubleshooting; reducing ``mixing_fraction``
makes either method more conservative.

After convergence, the Hamiltonian is solved once more using the converged
density. The returned orbitals, density, orbital energies, and energy
components are then recomputed as one final post-convergence result. They are
mutually consistent to the requested SCF tolerance rather than describing an
intermediate mixed density.

The reported total energy is assembled as

.. math::

   E_\mathrm{tot} = E_\mathrm{kin}
   + E_\mathrm{electron-ion}^{\mathrm{local}}
   + E_\mathrm{electron-ion}^{\mathrm{nonlocal}}
   + E_\mathrm{H}
   + E_\mathrm{xc}
   + E_\mathrm{ion-ion}.

Orbitals, phases, and degeneracy
--------------------------------

Plane waves and numerical eigenvectors are complex even at the
:math:`\Gamma` point. A non-degenerate time-reversal-symmetric orbital can
normally be made real by a single global phase rotation. Degenerate subspaces
are different: the eigensolver may return complex unitary mixtures of equally
valid real orbitals. A global phase applied to each state separately cannot, in
general, undo such mixing.

The plotting code removes the best global phase before displaying the real
component. It can report the residual imaginary norm, but suppressing that
diagnostic does not change the underlying orbital. Central molecular planes can
also be nodal planes. Since independently normalized contours can magnify tiny
numerical residues, ``plane="auto"`` is useful when plotting a heterogeneous
set of orbitals.

Current scope and limitations
-----------------------------

The implementation intentionally has a focused scope:

* cubic unit cells only;
* :math:`\Gamma`-point sampling only;
* closed-shell, spin-unpolarized systems with an even electron count;
* integer occupations without smearing;
* LDA/VWN5 and PBE exchange-correlation functionals;
* no geometry optimization, forces, stresses, molecular dynamics, or response
  properties; and
* a common working FFT grid rather than a production-style dual-grid scheme.

These constraints keep the numerical path compact enough to study while still
including spherical plane-wave cutoffs, iterative diagonalization, periodic
electrostatics, local and non-local pseudopotentials, and GPU execution.

Further reading
---------------

* R. M. Martin, *Electronic Structure: Basic Theory and Practical Methods*,
  Cambridge University Press.
* M. C. Payne *et al.*, "Iterative minimization techniques for ab initio
  total-energy calculations," *Reviews of Modern Physics* **64**, 1045 (1992).
* S. Goedecker, M. Teter, and J. Hutter, "Separable dual-space Gaussian
  pseudopotentials," *Physical Review B* **54**, 1703 (1996).
* J. P. Perdew, K. Burke, and M. Ernzerhof, "Generalized Gradient Approximation
  Made Simple," *Physical Review Letters* **77**, 3865 (1996).
* An open-access introduction to Hartree--Fock and localized-orbital DFT is
  available in `Elements of Electronic Structure Theory
  <https://ifilot.pages.tue.nl/elements-of-electronic-structure-theory/>`_.

Other educational software
--------------------------

* `PyQInt <https://pyqint.imc-tue.nl/>`_: Hartree--Fock and Gaussian-integral
  code.
* `PyDFT <https://pydft.imc-tue.nl/>`_: localized-orbital density functional
  theory code.
* `HFCXX <https://github.com/ifilot/hfcxx>`_: Hartree--Fock code in C++.
* `DFTCXX <https://github.com/ifilot/dftcxx>`_: localized-orbital DFT code in
  C++.
