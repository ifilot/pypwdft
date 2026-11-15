import unittest
import pytest
import numpy as np
import scipy as sc
import pyfftw

class TestFFT(unittest.TestCase):
    """
    Test whether NumPy FFT and PyFFTW give the same result
    """

    def test_fft(self):
        np.random.seed(42)
        sh = (16, 16, 16)
        data = np.random.random(sh)
        
        a = pyfftw.empty_aligned(sh, dtype='complex128')
        b = pyfftw.empty_aligned(sh, dtype='complex128')
        
        # Plan an fft over the last axis
        fft_object_a = pyfftw.builders.fftn(a)
        
        # perform FFT
        out = fft_object_a(data)
        
        # calculate FFT using numpy
        res = np.fft.fftn(data)
        
        np.testing.assert_almost_equal(out, res)
        
    def test_fft_conv(self):
        np.random.seed(42)
        sz  = 16
        sh = (sz,sz,sz)
        
        # create original data and its Fourier transform
        data = np.random.random(sh)
        fft_data = np.fft.ifftn(data)

        # produce convolution of FFT coefficients
        fft_conv = np.zeros_like(fft_data)
        sc.ndimage.convolve(fft_data, fft_data, fft_conv, mode='wrap')
        
        # cast back to realspace
        data2 = np.real(np.fft.fftn(fft_conv))

        # assert equality
        np.testing.assert_almost_equal(np.abs(data2), data**2)

if __name__ == '__main__':
    unittest.main()


@pytest.mark.parametrize("backend", ["pyfftw", "scipy", "numpy"])
@pytest.mark.parametrize("threads", [1, 2])
def test_grid_fft_matches_numpy(backend, threads):
    """
    Every CPU backend of GridFFT reproduces numpy.fft, for real and complex
    input, and leaves its input untouched.
    """
    from pypwdft.fft import GridFFT

    rng = np.random.default_rng(0)
    real = rng.random((6, 6, 6))
    field = real + 1j * rng.random((6, 6, 6))
    original = field.copy()
    fft = GridFFT(backend, threads=threads)

    np.testing.assert_allclose(fft.fftn(field), np.fft.fftn(field), atol=1e-12)
    np.testing.assert_allclose(fft.ifftn(field), np.fft.ifftn(field), atol=1e-12)
    np.testing.assert_allclose(fft.fftn(real), np.fft.fftn(real), atol=1e-12)
    np.testing.assert_array_equal(field, original)

    # results are independent copies, not views of a reused plan buffer
    first = fft.fftn(field)
    reference = first.copy()
    fft.fftn(real)
    np.testing.assert_array_equal(first, reference)


@pytest.mark.parametrize("threads", [0, -1, 1.5, True])
def test_invalid_thread_count(threads):
    from pypwdft.fft import GridFFT

    with pytest.raises(ValueError):
        GridFFT("pyfftw", threads=threads)


def test_openmp_limit_warning(monkeypatch):
    """A thread request above the OpenMP limit is reported, not ignored."""
    import pypwdft.fft as fft

    monkeypatch.setattr(fft.pyfftw, "_threading_type", "OMP", raising=False)
    monkeypatch.setattr(fft, "openmp_thread_limit", lambda: 1)
    with pytest.warns(RuntimeWarning, match="OpenMP is limited"):
        fft.check_fftw_threads(4)
