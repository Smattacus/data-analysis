#A few methods for getting normalized spectra.
#The scipy definitions for FFT involve no normalization forwards and 1 / N
#normalization backwards. These methods change it so that we normalize with
#1/sqrt(N) in each direction and also generates a frequency array alongside it.

import scipy as sp
from scipy import fftpack as fp

#These functions assume 1d arrays.

def spec(x, dt):
    """
Returns the normalized power spectrum (1 / sqrt(n) forward, 1 / sqrt(n)
backwards). This is so that Perseval's rule may be used easily.
    """
    n = size(x)
    f0 = 1/abs(n * dt)
    P = sp.fftpack.fft(x)
    f = (linspace(0, n, n) - n/2) * f0
    f = fp.fftshift(f)
    P = fp.fftshift(P)
    return f, P
    
    
    
