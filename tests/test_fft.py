import unittest
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
