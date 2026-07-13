import numpy as np
import pytest

from pypwdft.gth import GTHPseudopotential
from pypwdft.psystem import PeriodicSystem
from pypwdft.pypwdft import PyPWDFT
from pypwdft.pypwdft import LinOpH


def carbon_system(position=(5, 5, 5)):
    system = PeriodicSystem(10, ecut=4)
    system.add_atom(*position, 6)
    return system


def test_gth_parser_and_valence_electrons():
    system = carbon_system()
    system.add_atom(6.4, 5, 5, 1)
    gth = GTHPseudopotential(system)

    assert gth.symbols == ("C", "H")
    np.testing.assert_array_equal(gth.valence_charges, [4, 1])
    assert gth.nelec == 5
    assert gth["C"]["Zion"] == 4
    assert gth["C"]["rloc"] == pytest.approx(0.34883045)
    np.testing.assert_allclose(gth["C"]["cloc"][:2], [-8.51377110, 1.22843203])


def test_bundled_pbe_family():
    gth = GTHPseudopotential(carbon_system(), family="pbe")

    assert gth.family == "pbe"
    assert gth["C"]["rloc"] == pytest.approx(0.33847124)


def test_nonlocal_operator_is_hermitian():
    system = carbon_system()
    gth = GTHPseudopotential(system)
    rng = np.random.default_rng(1234)
    left = rng.normal(size=system.get_n_plane_waves()) + 1j * rng.normal(
        size=system.get_n_plane_waves()
    )
    right = rng.normal(size=system.get_n_plane_waves()) + 1j * rng.normal(
        size=system.get_n_plane_waves()
    )

    lhs = np.vdot(left, gth.apply_nonlocal(right))
    rhs = np.vdot(gth.apply_nonlocal(left), right)
    assert lhs == pytest.approx(rhs, abs=1e-12)
    assert gth.nprojectors == 1


def test_nonlocal_energy_matches_expectation_value():
    system = carbon_system()
    gth = GTHPseudopotential(system)
    rng = np.random.default_rng(42)
    orbitals = rng.normal(size=(system.get_n_plane_waves(), 2))
    orbitals = np.linalg.qr(orbitals)[0]

    expected = 2 * np.real(
        np.sum(orbitals.conj() * gth.apply_nonlocal(orbitals))
    )
    assert gth.nonlocal_energy(orbitals) == pytest.approx(expected)


def test_nonlocal_operator_is_used_by_hamiltonian():
    system = carbon_system()
    gth = GTHPseudopotential(system)
    rng = np.random.default_rng(7)
    coefficients = rng.normal(size=system.get_n_plane_waves()) + 1j * rng.normal(
        size=system.get_n_plane_waves()
    )
    npts = system.get_npts()
    hamiltonian = LinOpH(
        np.zeros((npts, npts, npts)),
        npts,
        system.get_pw_k2(),
        fft="numpy",
        pw_mask=system.get_pw_mask(),
        nonlocal_operator=gth.apply_nonlocal,
    )

    kinetic = 0.5 * system.get_pw_k2().ravel()[system.get_pw_mask().ravel()]
    expected = kinetic * coefficients + gth.apply_nonlocal(coefficients)
    np.testing.assert_allclose(hamiltonian @ coefficients, expected)


def test_geometry_change_requires_rebuild():
    system = carbon_system()
    gth = GTHPseudopotential(system)
    system.translate((1, 0, 0))

    with pytest.raises(RuntimeError, match="rebuild"):
        gth.local_potential()
    gth.rebuild()
    assert gth.local_potential().shape == (system.get_npts(),) * 3


def test_ewald_accepts_explicit_ionic_charges():
    system = PeriodicSystem(10, ecut=1)
    system.add_atom(4.3, 5, 5, 6)
    system.add_atom(5.7, 5, 5, 1)

    explicit = system.calculate_ewald_sum(charges=np.array([4, 1]))
    assert explicit != pytest.approx(system.calculate_ewald_sum())


def test_gth_scf_integration():
    system = PeriodicSystem(10, ecut=2)
    system.add_atom(4.3, 5, 5, 1)
    system.add_atom(5.7, 5, 5, 1)
    gth = GTHPseudopotential(system)

    result = PyPWDFT(
        system, fft="numpy", pseudopotential=gth
    ).scf(tol=1e-4, density_tol=1e-3)

    assert result["converged"]
    assert result["pseudopotential"] is gth
    assert result["Enonloc"] == 0
    assert result["Etot"] == pytest.approx(
        result["Ekin"]
        + result["Enuc"]
        + result["Enonloc"]
        + result["Erep"]
        + result["Eewald"]
        + result["Exc"]
    )


@pytest.mark.e2e
def test_pbe_gth_scf_integration():
    system = PeriodicSystem(10, ecut=2)
    system.add_atom(4.3, 5, 5, 1)
    system.add_atom(5.7, 5, 5, 1)
    gth = GTHPseudopotential(system, family="pbe")

    result = PyPWDFT(
        system,
        fft="numpy",
        functional="pbe",
        pseudopotential=gth,
    ).scf(tol=1e-6, density_tol=1e-4, maxiter=150)

    assert result["converged"]
    assert result["functional"] == "pbe"
    assert result["Etot"] == pytest.approx(-0.9687772867266564, abs=1e-8)
