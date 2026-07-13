import numpy as np
import pytest

from pypwdft.orbitals import imaginary_fraction, rotate_wavefunction_real


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
