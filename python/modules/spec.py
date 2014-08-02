#A few methods for getting normalized spectra.
#The scipy definitions for FFT involve no normalization forwards and 1 / N
#normalization backwards. These methods change it so that we normalize with
#1/sqrt(N) in each direction and also generates a frequency array alongside it.

import scipy as sp
from scipy import fftpack as fp
from numpy import size, linspace

#These functions assume 1d arrays.

def spec(x, dt):
    """
Returns the normalized power spectrum (1 / sqrt(n) forward, 1 / sqrt(n)
backwards). This is so that Perseval's rule may be used easily.
Usage:
    [f, g] = spec(x, dt)
INPUTS:
    x    - Array of data to perform the FFT on.
    dt   - Time element between array elements. (Scalar). Used to generate the
            frequency array.
OUTPUTS:
    f   - Frequency Array
    g   - FFT array.
    """
    n = size(x)
    f0 = 1/abs(n * dt)
    P = sp.fftpack.fft(x)
    #Generate an array 1:N, subtract (N + 1) /2 if odd, (N+2)/2 if even.
    f = (linspace(1, n, n) - (n + (1 + (n % 2 ==0 )))/2) * f0
    P = fp.fftshift(P)
    return f, P
    
    
    
