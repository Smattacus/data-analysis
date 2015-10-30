#This is a set of modules for simulating the data that is acquired in a typical 
#laser experiment. For now, basic modeling is done via Poisson distributions.

import spec
import numpy as np
import numpy.random as rand
import scipy.signal as sig

#For now, start by creating already  demodulated data.
def demoddata(bglam, flucf, fluca, t, seed=None):
    """
        Generates a demodulated data set with an underlying fluctuation.
        Generates background plasma noise according to the Poisson distribution
        with mean bglam, and a fluctuation signal with amplitude fluca and
        frequency flucf.
        
        INPUTS:
        bglam       :   Mean of poisson process of background counts.
        flucf       :   Frequency of fluctuation signal.
        fluca       :   Amplitude of fluctuation signal.

        Generates a sine wave that is added to a background Poisson process. By
        adjusting bglam and fluca, one may adjust the relative sizes of the
        background noise vs the sinusoidal signal.
    """
    bg = rand.poisson(lam = bglam, size=np.size(t))
    sig = fluca * np.sin(2 * np.pi * flucf * t)
    output = sig + bg
    return output

def uxcorr(x1, x2, *args):
    '''
        Uses fftconvolve with its unbiasing coefficient in order to create an
        unbiased cross correlation array. Assumes x1 and x2 are arrays.

        INPUTS:
        x1      :       First array to cross correlate.
        x2      :       Second array to cross correlate.
    '''
    #Reverse the last axis' ordering.
    xc = sig.fftconvolve(x1, x2[Ellipsis, ::-1], *args)
    N = (xc.shape[-1] + 1) / 2
    xc = xc / (N - np.abs(np.linspace(-(N-1), N-1, 2 * N - 1)))
    return xc

def averagexcorrs(navg, bg1, bg2, flucf, fluca1, fluca2, t):
    '''
        Performs an average of <navg> cross correlations between mean -
        subtracted PMT1 and PMT2. fluca1 and fluca2 determine the amplitude of
        the fluctuation in each PMT (this is a rough model for relative laser
        powers.)

        It is allowed that separate average background count rates are given.

        INPUTS:
        
        navg        :       Number of cross correlations to generate and
                            average.
        bg1         :       Mean for poisson count process for PMT 1.
        bg2         :       Mean for poisson count process for PMT2.
        flucf       :       Fluctuation fruquency.
        fluca1      :       Amplitude of fluctuation for pmt1.
        fluca2      :       Amplitude of fluctuation for pmt2.
        t           :       Time array for data. This models the demodulated
                            time axis.
        OUTPUTS:

        xcav        :       Averaged cross correlations.
    '''
    N = np.size(t)
    xcna = np.zeros((2 * N - 1, navg))
    for i in range(navg):
        pmt1 = demoddata(bg1, flucf, fluca1, t)
        pmt2 = demoddata(bg2, flucf, fluca2, t)
        pmt1 -= np.mean(pmt1)
        pmt2 -= np.mean(pmt2)
        xcna[:,i] = uxcorr(pmt1, pmt2)
    return np.mean(xcna, axis=1)

def axc_ndarr(navg, dims, bgs, flucfs1, flucfs2, flucas1, flucas2, t):
    '''
    Routine for generating a two dimensional array of averaged cross
    correlations.

    INPUTS:
    navg        :       Number of cross correlations to generate and average at
    each point.
    dims        :       Dimensions of the array to generate xcorrs for.
                        Averaged xcorrs will be generated for the last point
    bgs         :       Background count rate. Assumed constant for all
                        elements.
    flucfs1     :       Array of dimensions dims with fluctuation frequencies.
                        Each element is a tuple with (flucf1, flucf2).
    flucfs2     :       Array of dimensions dims with fluctuation frequencies
                        for PMT 2.
    flucas1     :       Array of dimensionality dims with fluctuation
                        amplitudes for PMT 1. Each element is a tuple with (fluca1,
                        fluca2)
    flucas2     :       Array of dimensionality dims with fluctuation
                        amplitudes for PMT 2.

    As an exmaple, if I wanted to generate a 2x2 matrix of cross correlations,
    each of which had been averaged over 50 runs measuring a constant 1200 Hz
    fluctuation, I'd say:

    bgs = 10000
    flucfs = np.array([1200,1200])
    dims = (2, 2)
    navg = 50
    xc2d = axc_2darr(50, dims, bgs, flucfs, flucas, t)
    '''
    #Generate all the background arrays.
    dims = dims + (navg,)
    bg1 = rand.poisson(bgs, size=dims + (np.size(t),))
    bg2 = rand.poisson(bgs, size=dims + (np.size(t),))
    #Generate the ndimensional PMT1 and PMT2 matrices.
    #Generate meshgrids for each of the desired coordinates (t, flucas1,
    #flucas2, flucfs1, flucfs2)
    #Each one of these meshgrid lines generates an ndimensional grid of what we
    #want to use for making the function.
    fa1g, fs1g, xcpref, tt = np.meshgrid(*[flucfs1, flucas1, np.ones((navg,)), t])
    fa2g, fs2g, xcpref, tt = np.meshgrid(*[flucfs2, flucas2, np.ones((navg,)), t])
    PMT1 = fa1g * np.sin(2 * np.pi * fs1g * tt) + bg1
    PMT2 = fa2g * np.sin(2 * np.pi * fs2g * tt) + bg2
    #Add the background noise to it.
    #Cross correlate the two along each last dimension.
    #Use a helper function which cross correlates the corresponding 1d arrays
    #in a certain dimension.
    uxc12 = uxcorr_ndim(PMT1, PMT2, axis=-1)
    return (uxc12, PMT1, PMT2)

def uxcorr_ndim(x1, x2, axis=-1):
    p12j = np.concatenate((x1, x2), axis=axis)
    return np.apply_along_axis(uxcorr_axis, axis, p12j)
    
def uxcorr_axis(arr):
    '''
    Helper function to be used with np.apply_along_axis.
    '''
    x1 = arr[0:arr.shape[0] / 2]
    x2 = arr[arr.shape[0]/2:]
    print(x1)
    print(x2)
    print(uxcorr(x1,x2))
    return uxcorr(x1, x2)

