#A few methods for getting normalized spectra.
#The scipy definitions for FFT involve no normalization forwards and 1 / N
#normalization backwards. These methods change it so that we normalize with
#1/sqrt(N) in each direction and also generates a frequency array alongside it.

import scipy as sp
from scipy import fftpack as fp
from numpy import size, linspace
from math import sqrt

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
    P = sp.fftpack.fft(x)/sqrt(n)
    #Generate an array 1:N, subtract (N + 1) /2 if odd, (N+2)/2 if even.
    f = (linspace(1, n, n) - (n + (1 + (n % 2 ==0 )))/2) * f0
    P = fp.fftshift(P)
    return f, P

def spec_dct(x, dt):
    """
Returns the normalized power spectrum (1 / sqrt(n) forward, 1 / sqrt(n)
backwards). This is so that Perseval's rule may be used easily. This particular
method uses the dct, which assumes the spectrum is real valued and even.
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
    P = sp.fftpack.rfft(x)/sqrt(n)
    #Generate an array 1:N, subtract (N + 1) /2 if odd, (N+2)/2 if even.
    f = (linspace(1, n, n) - (n + (1 + (n % 2 ==0 )))/2) * f0
    P = fp.fftshift(P)
    return f, P

def ispec(x, df):
    """
    Takes the inverse of a spectrum as calculated by spec().

    Usage:
    [t, g] = spec(x, dt)
INPUTS:
    x   - array of data to inverse transform
    df  - Frequency elements between array elements. (Scalar). Used to
    generate the frequency array.
OUTPUTS:
    t   - Time array.
    g   - Real - time function.
"""
    n = size(x)
    T = 1/abs(df)
    g = sp.fftpack.ifft(fp.ifftshift(x)) * sqrt(n)
    t = (linspace(0, n-1, n)) * T / n
    return t, g
    
    
#2D arrays.

#2D version of spec,
def spec2d(x, dtx, dty):
    """Returns the normalized, fftshifted power spectrum. The spectrum is
    normalized according to 1/sqrt(x.shape[0]) 1 / sqrt(x.shape[1]). fft is
    performed along each row and column, and then each row and column is
    shifted so that the center of the array is the 0 frequency point.

    USAGE:
        [f1, f2, g] = spec(x, dt)

    INPUTS:
        x   - 2D numpy array of data to perform the fft2 on.
        dtx     - time element between x array elements. (Scalar).
        dty     - Time element betweern y array elements.
    OUTPUTS:
        f1   - Frequency array for axis 1.
        f2  - Frequency array for axis 2.
        g   - 2D FFT array.
    """
    n1 = x.shape[0]
    n2 = x.shape[1]
    f0_1 = 1 / abs(n1 * dtx)
    f0_2 = 1 / abs(n2 * dty)
    P = fp.fft2(x) / sqrt(n1) / sqrt(n2)
    f1 = (linspace(1, n1, n1) - (n1 + (1 + (n1 % 2 ==0 )))/2) * f0_1
    f2 = (linspace(1, n2, n2) - (n2 + (1 + (n2 % 2 ==0 )))/2) * f0_2
    P = fp.fftshift(P)
    return f1, f2, P

def ispec2d(g, df1, df2):
    """Returns the original, inverse FFTshifted power spectrum.

    USAGE:
        [t1, t2, x] = ispec2d(g, df1, df2)

    INPUTS:
        g   -   Spectra to invert.
        df1 -   Frequency step of first dimension.
        df2 -   Frequency step of second dimension.

    """
    n1 = g.shape[0]
    n2 = g.shape[1]
    X = fp.ifftshift(g)
    X = fp.ifft2(X) * sqrt(n1) * sqrt(n2)
    t0_1 = 1 / abs(n1 * df1)
    t0_2 = 1 / abs(n2 * df2)
    t1 = (linspace(1, n1, n1) - (n1 + (1 + (n1 % 2 ==0 )))/2) * t0_1
    t2 = (linspace(1, n1, n1) - (n1 + (1 + (n1 % 2 ==0 )))/2) * t0_2
    return t1, t2, X

