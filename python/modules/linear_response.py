#Some methods for calculating linear response functions
import numpy as np
from scipy import signal
from scipy import fftpack as fp
import spec
import Gnuplot
from scipy import signal


def response(xcm, a1, a2, tc, tw, dt):
    """Calculates and returns the linear response function corresponding to Cr
    and Cl, the folded cross correlation functions.

    (Hxy, Hyx) = response(xcm, a1, a2, tc, tw, dt)

    This method assumes that xcm was created with the matlab command xcorr as
    follows:
    xcm = xcorr(a1, a2, 'unbisaed');
    This results in Cr corresponding to the linear response of X to Y and Cl
    corresponding to the linear response of Y to X.

    This method uses a gaussian windowing function in order to clean up the
    signal before taking spectra. tc is a parameter to allow tuning of this window.

    Parameters:
    xcm : array_like
        Cross correlation of a1 and a2.
    a1 : array_like
        Autocorrelation of time series one.
    a2 : array_like
        Autocorelation of times series two.
    tc : float64
        Window width to use. Goes like win = exp(-(t/tc)**2/2)
    tw : float64
        Maximum time to include in data. All time points larger will be cut
        out. A typical value for a tc = 0.01 is 0.04.
    Returns:
    Hxy : array_like
        Linear response of X to Y.
    Hyx: array_like
        Linear response of Y to X.
    """
    N = (np.size(xcm) + 1) / 2
    CR = np.append(xcm[N-1:][::-1], xcm[N:])
    CL = np.append(xcm[0:N-1], xcm[0:N][::-1])
    t = np.linspace(-(N-1),(N-1), 2 * N - 1) / 1e5
    win = np.exp(-(t/tc)**2/2)
    #Window functions
    CRw = CR * win
    CLw = CL * win
    a1w = a1 * win
    a2w = a2 * win
    #Get the middle points according to tw
    tc = np.where(abs(t) < tw)
    CRw = CRw[tc]
    CLw = CLw[tc]
    a1w = a1w[tc]
    a2w = a2w[tc]
    #fftshift + reverse so that t = 0 is at the beginning (rather than the end)
    CRw = fp.fftshift(CRw)[::-1]
    CLw = fp.fftshift(CLw)[::-1]
    a1w = fp.fftshift(a1w)[::-1]
    a2w = fp.fftshift(a2w)[::-1]
    #take spectra
    [f, gCRw] = spec.spec(CRw, dt)
    [f, gCLw] = spec.spec(CLw, dt)
    [f, ga1w] = spec.spec(a1w, dt)
    [f, ga2w] = spec.spec(a2w, dt)
    #Get k space response functions
    #Really, using the reals here is kind of redundant since CR, CL, a1w, and
    #a2w are all symmetric by construction. Shifting it to t = 0 and taking
    #real() is just another way of getting abs(), only by making sure the
    #fft vector is on the real line. w/e
    gHxy = np.real(gCRw) / np.real(ga2w)
    gHyx = np.real(gCLw) / np.real(ga1w)
    gHxyHT = signal.hilbert(gHxy)
    gHyxHT = signal.hilbert(gHyx)
    [t, Hxy] = spec.ispec(gHxyHT, f[1] - f[0])
    [t, Hyx] = spec.ispec(gHyxHT, f[1] - f[0])
    return (Hxy[::-1], Hyx[::-1], t)

