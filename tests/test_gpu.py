import importlib.util
import sys

import numpy as np
import pytest

from pypwdft.gth import GTHPseudopotential
from pypwdft.psystem import PeriodicSystem
from pypwdft.pypwdft import PyPWDFT
from pypwdft.pypwdft import _load_cupy


def hydrogen_system():
    system = PeriodicSystem(10, ecut=2)
    system.add_atom(4.3, 5, 5, 1)
    system.add_atom(5.7, 5, 5, 1)
    return system


def test_unknown_fft_backend_is_rejected():
    with pytest.raises(ValueError, match="Unknown FFT backend"):
        PyPWDFT(hydrogen_system(), fft="not-a-backend")


def test_missing_cupy_has_actionable_error(monkeypatch):
    monkeypatch.setitem(sys.modules, "cupy", None)
    with pytest.raises(ImportError, match="requires a CuPy package"):
        _load_cupy()


def test_numpy_nonlocal_operator_factory_matches_original():
    system = PeriodicSystem(10, ecut=2)
    system.add_atom(5, 5, 5, 6)
    pseudopotential = GTHPseudopotential(system)
    rng = np.random.default_rng(91)
    coefficients = rng.normal(size=system.get_n_plane_waves())

    np.testing.assert_allclose(
        pseudopotential.nonlocal_operator(np)(coefficients),
        pseudopotential.apply_nonlocal(coefficients),
    )


def require_cuda():
    """Skip the calling test unless CuPy and a CUDA GPU are available."""
    if importlib.util.find_spec("cupy") is None:
        pytest.skip("CuPy is not installed")
    import cupy as cp

    try:
        if cp.cuda.runtime.getDeviceCount() < 1:
            pytest.skip("No CUDA GPU is available")
    except cp.cuda.runtime.CUDARuntimeError:
        pytest.skip("CUDA cannot be initialized")


def carbon_monoxide_scf(backend, eigensolver):
    """GTH-PBE calculation covering projectors, GGA and unoccupied bands."""
    system = PeriodicSystem(8, ecut=4)
    system.add_atom(2.934, 4, 4, 6)
    system.add_atom(5.066, 4, 4, 8)
    return PyPWDFT(
        system,
        fft=backend,
        pseudopotential=GTHPseudopotential(system, family="pbe"),
        functional="pbe",
    ).scf(tol=1e-7, nsol=6, eigensolver=eigensolver)


def test_eigensolvers_agree_on_cpu():
    davidson = carbon_monoxide_scf("numpy", "davidson")
    arpack = carbon_monoxide_scf("numpy", "arpack")

    assert davidson["Etot"] == pytest.approx(arpack["Etot"], abs=1e-8)
    np.testing.assert_allclose(davidson["orbe"], arpack["orbe"], atol=1e-5)
    assert (
        davidson["hamiltonian_applications"]
        < arpack["hamiltonian_applications"]
    )


@pytest.mark.gpu
@pytest.mark.parametrize("eigensolver", ["davidson", "arpack"])
def test_cupy_gth_pbe_scf_matches_numpy(eigensolver):
    require_cuda()

    cpu = carbon_monoxide_scf("numpy", eigensolver)
    gpu = carbon_monoxide_scf("cupy", eigensolver)

    assert gpu["Etot"] == pytest.approx(cpu["Etot"], abs=1e-8)
    np.testing.assert_allclose(gpu["orbe"], cpu["orbe"], atol=1e-5)
    np.testing.assert_allclose(gpu["edens"], cpu["edens"], atol=1e-5)


@pytest.mark.gpu
def test_cupy_scf_matches_numpy():
    require_cuda()

    cpu = PyPWDFT(hydrogen_system(), fft="numpy").scf(
        tol=1e-4, density_tol=1e-3
    )
    gpu = PyPWDFT(hydrogen_system(), fft="cupy").scf(
        tol=1e-4, density_tol=1e-3
    )

    assert gpu["fft"] == "cupy"
    assert isinstance(gpu["edens"], np.ndarray)
    assert isinstance(gpu["orbc_rs"], np.ndarray)
    assert gpu["Etot"] == pytest.approx(cpu["Etot"], abs=1e-7)
    np.testing.assert_allclose(gpu["orbe"], cpu["orbe"], atol=1e-7)
