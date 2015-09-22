#Functions for performing Wiener deconvolution.

import numpy as np

def FilterWhiteNoise(f, H, S, fwn):
    """
        Filters, using a wiener filter, a spectrum assuming that white noise is
        the primary source of noise. Does windowing to estimate the noise.

        INPUTS:
        f   -   Frequency of the spectra.
        S   -   Observed signal spectrum.
        H   -   Spectrum of the known signal which one wants to deconvolve.
        fwn -   Tuple of the ends of a range of frequencies to use for noise estimation (via mean). White
                    noise assumption. (fstart, fend)
    """
    fn = np.where((f < fwn[1]) - (f < fwn[0]))
    noise = np.mean(S[fn])
    G = WienerFilter(H, S, noise)
    DC = G * S
    return(DC)



def WienerFilter(H, S, N):
    """

    G = WienerFilter(H, S, N)

    Returns a Wiener Filter based on the spectra sent in to it. This filter
    can then be pointwise multiplied with the relevant spectrum (S) in order to 
    filter the resulting data set.

    INPUTS:
        H   -   Spectrum of the known part which one wants to deconvolve.
        S   -   Observed spectrum.
        N   -   Noise spectrum.
    OUTPUTS:
        G   -   Wiener deconvolution filter.
    """
    G = 1 / H * (np.abs(H)**2 / (np.abs(H)**2 + N / S))
    return G

def WienerFilter2D(xg1, xg2, S, H, N):
    '''
    Implements a wiener filter on a two dimensional matrix (generally one
    created with meshgrid)
    '''
