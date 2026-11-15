import numpy as np
import pytest

from pypwdft.davidson import conjugation_partner, davidson, initial_guess


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


def test_hamiltonian_block_is_split_into_bounded_batches(monkeypatch):
    """A block larger than the FFT batch limit gives the same result."""
    import pypwdft.pypwdft as engine
    from pypwdft.psystem import PeriodicSystem

    system = PeriodicSystem(8, ecut=3)
    system.add_atom(4, 4, 4, 6)
    npts = system.get_npts()
    rng = np.random.default_rng(7)
    potential = rng.normal(size=(npts, npts, npts))
    npw = system.get_n_plane_waves()
    block = rng.normal(size=(npw, 7)) + 1j * rng.normal(size=(npw, 7))

    def apply():
        return engine.LinOpH(
            potential, npts, system.get_pw_k2(), fft="numpy",
            pw_mask=system.get_pw_mask(),
        ).matmat(block)

    reference = apply()
    # room for four grids per batch: seven states become batches of 4, 2, 1
    monkeypatch.setattr(engine, "FFT_BATCH_BYTES", 4 * 16 * npts**3)
    np.testing.assert_allclose(apply(), reference, atol=1e-12)
    monkeypatch.setattr(engine, "FFT_BATCH_BYTES", 1)
    np.testing.assert_allclose(apply(), reference, atol=1e-12)


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


def real_space_real_problem(ecut=6, nonlocal_operator=False):
    """A Hamiltonian that is real in real space, with degenerate levels."""
    from pypwdft.gth import GTHPseudopotential
    from pypwdft.psystem import PeriodicSystem
    from pypwdft.pypwdft import LinOpH

    system = PeriodicSystem(10, ecut=ecut)
    npts = system.get_npts()
    axis = np.linspace(0, 10, npts, endpoint=False)
    # Spherical well on a grid point away from the origin and the cell
    # centre: exactly cubic, so its p-like level is exactly degenerate, yet
    # its eigenfunctions carry non-trivial plane-wave phases.
    centre = axis[[3, 5, 8]]
    system.add_atom(*centre, 6)
    z, y, x = np.meshgrid(axis, axis, axis, indexing="ij")
    distance2 = sum(
        ((coordinate - origin + 5) % 10 - 5)**2
        for coordinate, origin in zip((x, y, z), centre)
    )
    potential = -3 * np.exp(-distance2 / 4)
    mask = system.get_pw_mask()
    hamiltonian = LinOpH(
        potential, npts, system.get_pw_k2(), fft="numpy", pw_mask=mask,
        nonlocal_operator=(
            GTHPseudopotential(system).nonlocal_operator()
            if nonlocal_operator else None
        ),
    )
    return hamiltonian, 0.5 * system.get_pw_k2()[mask], mask


def test_conjugation_partner_maps_every_plane_wave_to_its_negative():
    from pypwdft.psystem import PeriodicSystem

    system = PeriodicSystem(10, ecut=6)
    mask = system.get_pw_mask()
    partner = conjugation_partner(mask)
    kvec = system.get_pw_k()[mask]

    np.testing.assert_allclose(kvec[partner], -kvec, atol=1e-12)
    np.testing.assert_array_equal(partner[partner], np.arange(len(partner)))
    # a mask that is not closed under G -> -G has no partner map
    lopsided = np.zeros((4, 4, 4), dtype=bool)
    lopsided[1, 0, 0] = True
    assert conjugation_partner(lopsided) is None


@pytest.mark.parametrize("nonlocal_operator", [False, True])
def test_hamiltonian_maps_real_functions_onto_real_functions(nonlocal_operator):
    """The premise of the real-orbital Davidson solver, non-local part included."""
    hamiltonian, _, mask = real_space_real_problem(
        nonlocal_operator=nonlocal_operator
    )
    partner = conjugation_partner(mask)
    rng = np.random.default_rng(1)
    block = rng.normal(size=(len(partner), 3)) + 1j * rng.normal(size=(len(partner), 3))
    block = 0.5 * (block + block[partner].conj())

    applied = hamiltonian.matmat(block)

    np.testing.assert_allclose(applied, applied[partner].conj(), atol=1e-12)


@pytest.mark.parametrize("nstates", [1, 2, 3, 4, 5, 7])
def test_initial_guess_of_real_functions_has_full_rank(nstates):
    _, kinetic, mask = real_space_real_problem()
    partner = conjugation_partner(mask)

    guess = initial_guess(kinetic, nstates, partner=partner)

    np.testing.assert_allclose(guess, guess[partner].conj(), atol=1e-15)
    assert np.linalg.matrix_rank(guess) == nstates


@pytest.mark.parametrize("nstates", [2, 3, 4, 5])
def test_davidson_returns_real_orbitals_also_within_a_cut_level(nstates):
    """
    Levels 2-4 are degenerate. Every eigenvector is a real function in real
    space, also when only part of the degenerate level is requested, and no
    phase rotation is needed.
    """
    hamiltonian, kinetic, mask = real_space_real_problem()
    partner = conjugation_partner(mask)

    energies, vectors, converged = davidson(
        hamiltonian.matmat, initial_guess(kinetic, nstates, partner=partner),
        kinetic, tol=1e-9, partner=partner,
    )
    reference, _, _ = davidson(
        hamiltonian.matmat, initial_guess(kinetic, 5), kinetic, tol=1e-9
    )

    assert converged
    np.testing.assert_allclose(energies, reference[:nstates], atol=1e-9)
    assert reference[3] - reference[1] < 1e-8          # the degenerate level
    np.testing.assert_allclose(vectors, vectors[partner].conj(), atol=1e-12)
    grid = np.zeros((vectors.shape[1], mask.size), dtype=complex)
    grid[:, mask.ravel()] = vectors.T
    orbitals = np.fft.ifftn(grid.reshape((-1,) + mask.shape), axes=(1, 2, 3))
    assert np.abs(orbitals.imag).max() < 1e-12 * np.abs(orbitals).max()
    np.testing.assert_allclose(
        vectors.conj().T @ vectors, np.eye(nstates), atol=1e-12
    )
