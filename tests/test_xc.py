import numpy as np
import pytest

from pypwdft import PeriodicSystem, PyPWDFT


def test_vwn5_correlation_reference_values():
    """Check the spin-unpolarized VWN5 Ceperley-Alder fit."""
    calculator = PyPWDFT(PeriodicSystem(10, ecut=1), fft="numpy")
    rs = np.array([0.5, 1.0, 2.0, 5.0, 10.0])
    density = 3 / (4 * np.pi * rs**3)

    energy, potential = calculator._PyPWDFT__lda_c_vwn(density)

    np.testing.assert_allclose(
        energy,
        [
            -0.07706330702344717,
            -0.060018686442541096,
            -0.044782788614621816,
            -0.028133762289731405,
            -0.01854452716940295,
        ],
        rtol=1e-13,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        potential,
        [
            -0.08562449002101037,
            -0.06781621037986249,
            -0.05160382394979036,
            -0.03338417103536458,
            -0.022518326145863098,
        ],
        rtol=1e-13,
        atol=1e-15,
    )


def test_pbe_periodic_reference_values():
    """Check PBE energy and its FFT-based GGA potential."""
    system = PeriodicSystem(10, ecut=1)
    calculator = PyPWDFT(system, fft="numpy", functional="pbe")
    points = system.get_r()
    density = (
        0.02
        + 0.003 * np.cos(2 * np.pi * points[..., 0] / 10)
        + 0.002 * np.sin(4 * np.pi * points[..., 1] / 10)
    )

    energy, potential = calculator._PyPWDFT__calculate_xc(density)

    np.testing.assert_allclose(
        energy.ravel()[:5],
        [
            -0.25305940197816534,
            -0.25069659843863960,
            -0.24447546464835904,
            -0.23693137042042411,
            -0.23169446455810155,
        ],
        rtol=1e-13,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        potential.ravel()[:5],
        [
            -0.32968259704856270,
            -0.32658838210162240,
            -0.31843547726467614,
            -0.30853389968516737,
            -0.30166431680805283,
        ],
        rtol=1e-13,
        atol=1e-15,
    )


def test_unknown_functional_is_rejected():
    with pytest.raises(ValueError, match="Unknown exchange-correlation"):
        PyPWDFT(PeriodicSystem(10, ecut=1), functional="not-a-functional")
