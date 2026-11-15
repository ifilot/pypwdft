import numpy as np
import pytest

from pypwdft.orbitals import (
    degenerate_groups,
    imaginary_fraction,
    rotate_orbitals_real,
    rotate_wavefunction_real,
)


def test_rotate_wavefunction_removes_a_global_phase():
    wavefunction = np.array([1.0, -2.0, 3.0]) * np.exp(0.73j)

    rotated = rotate_wavefunction_real(wavefunction)

    assert rotated.real == pytest.approx([1.0, -2.0, 3.0])
    assert rotated.imag == pytest.approx(0.0, abs=1e-14)
    assert np.abs(rotated) == pytest.approx(np.abs(wavefunction))


def test_rotate_wavefunction_minimizes_imaginary_norm():
    wavefunction = np.array([1.0 + 2.0j, -0.5 + 0.25j, 3.0 - 1.0j])
    rotated = rotate_wavefunction_real(wavefunction)
    trial_norms = [
        np.linalg.norm((wavefunction * np.exp(-1j * phase)).imag)
        for phase in np.linspace(0.0, np.pi, 2001)
    ]

    assert np.linalg.norm(rotated.imag) <= min(trial_norms) + 1e-6


def test_imaginary_fraction_handles_real_and_zero_fields():
    assert imaginary_fraction(np.array([1.0, 2.0])) == 0.0
    assert imaginary_fraction(np.zeros(3, dtype=complex)) == 0.0


def real_orthonormal_set(count, size=400, seed=0):
    """Real orthonormal functions, stored as rows."""
    rng = np.random.default_rng(seed)
    return np.linalg.qr(rng.normal(size=(size, count)))[0].T


def random_unitary(size, seed=0):
    rng = np.random.default_rng(seed)
    matrix = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    return np.linalg.qr(matrix)[0]


def test_degenerate_groups_follow_energy_gaps():
    groups = degenerate_groups([-1.0, -0.5, -0.5 + 1e-7, 0.2, 0.2, 0.2, 0.9])

    assert [list(group) for group in groups] == [[0], [1, 2], [3, 4, 5], [6]]
    assert degenerate_groups([]) == []


@pytest.mark.parametrize("degeneracy", [2, 3, 5])
def test_complex_mixture_of_a_degenerate_level_becomes_real(degeneracy):
    """
    A unitary mixture such as (p_x + i p_y)/sqrt(2) cannot be made real by a
    phase per orbital, but the level as a whole can.
    """
    real = real_orthonormal_set(degeneracy + 2)
    energies = np.array([-1.0] + [-0.3] * degeneracy + [0.4])
    mixed = real.astype(complex)
    mixed[1:-1] = random_unitary(degeneracy, seed=degeneracy) @ real[1:-1]
    mixed[0] *= np.exp(0.4j)
    mixed[-1] *= np.exp(-1.9j)

    per_orbital = [
        imaginary_fraction(rotate_wavefunction_real(orbital))
        for orbital in mixed[1:-1]
    ]
    assert max(per_orbital) > 0.05

    rotated = rotate_orbitals_real(mixed, energies)

    assert max(imaginary_fraction(orbital) for orbital in rotated) < 1e-12
    # still orthonormal, and every level still spans the same space
    assert rotated.conj() @ rotated.T == pytest.approx(
        np.eye(len(rotated)), abs=1e-12
    )
    overlap = rotated[1:-1].conj() @ real[1:-1].T
    assert np.linalg.svd(overlap, compute_uv=False) == pytest.approx(1.0)
    assert np.abs(rotated[0]) == pytest.approx(np.abs(real[0]))


def test_orbitals_that_are_real_up_to_a_phase_are_not_mixed():
    """Members of a level that only carry a phase keep their identity."""
    real = real_orthonormal_set(3, seed=4)
    phases = np.exp(1j * np.array([0.3, -1.1, 2.5]))
    rotated = rotate_orbitals_real(phases[:, None] * real, [0.1, 0.1, 0.1])

    assert rotated.imag == pytest.approx(0.0, abs=1e-12)
    assert np.abs(rotated.real) == pytest.approx(np.abs(real))


def test_single_orbitals_match_the_global_phase_rotation():
    rng = np.random.default_rng(2)
    orbitals = rng.normal(size=(3, 4, 4, 4)) + 1j * rng.normal(size=(3, 4, 4, 4))

    rotated = rotate_orbitals_real(orbitals, [0.0, 1.0, 2.0])
    without_energies = rotate_orbitals_real(orbitals)

    for index, orbital in enumerate(orbitals):
        expected = rotate_wavefunction_real(orbital)
        assert rotated[index] == pytest.approx(expected)
        assert without_energies[index] == pytest.approx(expected)
    assert rotated.shape == orbitals.shape


def test_incomplete_level_is_rotated_as_far_as_possible():
    """One member of a complex pair cannot be made real, and is reported."""
    real = real_orthonormal_set(2, seed=6)
    lone = (real[0] + 1j * real[1]) / np.sqrt(2)

    rotated = rotate_orbitals_real(lone[None], [0.0])

    assert np.linalg.norm(rotated) == pytest.approx(1.0)
    assert imaginary_fraction(rotated[0]) == pytest.approx(np.sqrt(0.5))


def test_rotate_orbitals_validates_energies():
    with pytest.raises(ValueError):
        rotate_orbitals_real(np.zeros((2, 3)), [0.0])
