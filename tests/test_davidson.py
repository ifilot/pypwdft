import numpy as np
import pytest

from pypwdft.davidson import davidson, initial_guess


def hermitian_problem(size=120, seed=3):
    """A dense Hermitian matrix dominated by a kinetic-like diagonal."""
    rng = np.random.default_rng(seed)
    kinetic = np.sort(rng.uniform(0, 20, size))
    coupling = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    matrix = np.diag(kinetic) + 0.05 * (coupling + coupling.conj().T)
    return matrix, kinetic


def test_davidson_matches_dense_diagonalisation():
    matrix, kinetic = hermitian_problem()
    nstates = 5

    energies, vectors, converged = davidson(
        lambda v: matrix @ v, initial_guess(kinetic, nstates), kinetic, tol=1e-9
    )

    assert converged
    np.testing.assert_allclose(
        energies, np.linalg.eigvalsh(matrix)[:nstates], atol=1e-12
    )
    np.testing.assert_allclose(
        vectors.conj().T @ vectors, np.eye(nstates), atol=1e-12
    )
    residual = matrix @ vectors - vectors * energies
    assert np.linalg.norm(residual, axis=0).max() < 1e-9


def test_hamiltonian_block_application_matches_single_vectors():
    from pypwdft.gth import GTHPseudopotential
    from pypwdft.psystem import PeriodicSystem
    from pypwdft.pypwdft import LinOpH

    system = PeriodicSystem(8, ecut=3)
    system.add_atom(4, 4, 4, 6)
    npts = system.get_npts()
    rng = np.random.default_rng(5)
    potential = rng.normal(size=(npts, npts, npts))
    npw = system.get_n_plane_waves()
    block = rng.normal(size=(npw, 3)) + 1j * rng.normal(size=(npw, 3))

    results = []
    for backend in ("numpy", "scipy", "pyfftw"):
        hamiltonian = LinOpH(
            potential, npts, system.get_pw_k2(), fft=backend,
            pw_mask=system.get_pw_mask(),
            nonlocal_operator=GTHPseudopotential(system).nonlocal_operator(),
        )
        single = np.stack(
            [hamiltonian.matvec(block[:, i]) for i in range(3)], axis=1
        )
        np.testing.assert_allclose(hamiltonian.matmat(block), single, atol=1e-12)
        assert hamiltonian.applications == 6
        results.append(single)
    np.testing.assert_allclose(results[0], results[1], atol=1e-12)
    np.testing.assert_allclose(results[0], results[2], atol=1e-12)


def test_davidson_reports_missing_convergence():
    matrix, kinetic = hermitian_problem()

    *_, converged = davidson(
        lambda v: matrix @ v, initial_guess(kinetic, 3), kinetic,
        tol=1e-12, maxiter=1,
    )

    assert not converged


def test_unknown_eigensolver_is_rejected():
    from pypwdft import SCFSettings

    with pytest.raises(ValueError, match="eigensolver"):
        SCFSettings(eigensolver="lanczos")
