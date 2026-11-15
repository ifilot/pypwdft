import numpy as np
import pytest
import scipy.sparse.linalg

from pypwdft.lanczos import lanczos_lowest, start_vector
from pypwdft.psystem import PeriodicSystem
from pypwdft.pypwdft import LinOpH


@pytest.fixture(scope="module")
def degenerate_hamiltonian():
    """
    A particle in a spherical well at the centre of a cubic cell. Its second
    level is triply degenerate (p-like), and its states have well-defined
    symmetries: the situation in which a plain Lanczos solve fails.
    """
    system = PeriodicSystem(10, ecut=10)
    system.add_atom(5, 5, 5, 6)
    npts = system.get_npts()
    axis = np.linspace(0, 10, npts, endpoint=False)
    z, y, x = np.meshgrid(axis, axis, axis, indexing="ij")
    potential = -3 * np.exp(-((x - 5)**2 + (y - 5)**2 + (z - 5)**2) / 4)
    mask = system.get_pw_mask()
    k2 = system.get_pw_k2()
    hamiltonian = LinOpH(potential, npts, k2, fft="numpy", pw_mask=mask)
    return hamiltonian, 0.5 * k2[mask]


LEVELS = [-1.408796, -0.544391, -0.544391, -0.544391, -0.120850]


def test_plain_lanczos_misses_degenerate_states(degenerate_hamiltonian):
    """Documents the failure that lanczos_lowest guards against."""
    hamiltonian, kinetic = degenerate_hamiltonian
    operator = scipy.sparse.linalg.LinearOperator(
        hamiltonian.shape, matvec=hamiltonian.matvec, dtype=complex
    )
    energies, _ = scipy.sparse.linalg.eigsh(
        operator, 4, which="SA", v0=start_vector(kinetic), tol=1e-6
    )

    assert not np.allclose(np.sort(energies), LEVELS[:4], atol=1e-3)


@pytest.mark.parametrize("nstates", [1, 2, 3, 4, 5])
@pytest.mark.parametrize("tol", [1e-2, 1e-6, 0])
def test_lanczos_lowest_finds_every_degenerate_state(
        degenerate_hamiltonian, nstates, tol):
    """
    Correct for every tolerance and for every number of states, including
    those that cut through the triply degenerate level.
    """
    hamiltonian, kinetic = degenerate_hamiltonian

    energies, vectors = lanczos_lowest(
        scipy.sparse.linalg, hamiltonian.matmat, nstates,
        start_vector(kinetic), tol=tol,
    )

    np.testing.assert_allclose(energies, LEVELS[:nstates], atol=1e-5)
    np.testing.assert_allclose(
        vectors.conj().T @ vectors, np.eye(nstates), atol=1e-10
    )
    residual = hamiltonian.matmat(vectors) - vectors * energies
    assert np.linalg.norm(residual, axis=0).max() < 1e-3


def test_start_vector_is_deterministic_and_generic():
    kinetic = np.linspace(0, 20, 50)

    vector = start_vector(kinetic)

    np.testing.assert_array_equal(vector, start_vector(kinetic))
    assert np.linalg.norm(vector) == pytest.approx(1.0)
    assert np.all(np.abs(vector) > 0)


def test_lanczos_reports_missing_convergence(monkeypatch):
    """A solve that never converges ends in an error, not in a long wait."""
    matrix = np.diag(np.arange(40.0))

    def never_converges(*args, **kwargs):
        raise scipy.sparse.linalg.ArpackNoConvergence("no convergence", None, None)

    monkeypatch.setattr(scipy.sparse.linalg, "eigsh", never_converges)
    with pytest.raises(RuntimeError, match="did not converge"):
        lanczos_lowest(
            scipy.sparse.linalg, lambda block: matrix @ block, 3,
            start_vector(np.arange(40.0)),
        )


@pytest.mark.parametrize("nstates", [1, 2, 3, 4, 5])
def test_lanczos_lowest_returns_real_orbitals(degenerate_hamiltonian, nstates):
    """
    With the -G partner map the eigenvectors are real functions in real
    space, also when ``nstates`` cuts through the degenerate level.
    """
    from pypwdft.davidson import conjugation_partner

    hamiltonian, kinetic = degenerate_hamiltonian
    mask = np.zeros(hamiltonian.npts**3, dtype=bool)
    mask[np.asarray(hamiltonian.pw_mask)] = True
    partner = conjugation_partner(mask.reshape((hamiltonian.npts,) * 3))

    energies, vectors = lanczos_lowest(
        scipy.sparse.linalg, hamiltonian.matmat, nstates,
        start_vector(kinetic), tol=1e-8, partner=partner,
    )

    np.testing.assert_allclose(energies, LEVELS[:nstates], atol=1e-5)
    np.testing.assert_allclose(vectors, vectors[partner].conj(), atol=1e-12)
    np.testing.assert_allclose(
        vectors.conj().T @ vectors, np.eye(nstates), atol=1e-10
    )
