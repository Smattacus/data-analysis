#This is a file containing major functions from MATLAB used for the purpose of
#reducing data from the Skiff lab HLMX experiment.

#This is to be used to replace the MATLAB stuff finally.

import os
import numpy as np
import scipy as sp
from scipy import signal as sig
from scipy import fftpack as fp
import glob
import h5py
import spec


def genAllTPBoxFindPhase(path, Fc, Fa):
    '''
Function to generate ALL the downsampled TOP and BOT arrays of LIF + Noise
and noise data. 

[TOP1, BOT1, TOP2, BOT2] = genAllTPBoxFindPhase(path, Fc, Fa)

This function assumes that the same phase is to be used across all files
in the directory pointed to by path. It generates a straight square wave
with phase and frequency given by phase and freq, respectively. This
square wave is used to demodulate the data.

This function finds all available files by searching for .h5 files. It
will use all of them in the directory, so delete or move unwanted ones.

INPUTS:
   path - Path to a directory of h5 files.
   Fc - Frequency of the square wave function. Generally, this is the
            laser chop frequency
   Fa   - Acquisition frequency.

OUTPUTS:
   TOP1 - NxM array of LIF +plasma data. N = number of files, M = # of points, PMT 1.
   BOT1 - NxM array of plasma data, PMT 1.
   TOP2 - NxM array of LIF + plasma data, PMT 2.
   BOT2 - NxM array of plasma data, PMT 2.
    '''
    os.chdir(path)
    l = glob.glob('*.h5')
    l.sort()
    
    nfiles = np.size(l)
    #'r' flag opens as read only. 'w' will overwrite the file!!
    f = h5py.File(l[0], 'r')
    d = np.array(f['PMT_DATA_8BIT'])
    #Python is row x column (row major), so the index for np.size is different than MATLAB's
    N = np.size(d, 0)
    total_t = N / Fa
    base_phase = np.linspace(0 , 2 * np.pi * total_t * Fc, N)
    TOP1 = np.zeros((nfiles, N * Fc / Fa))
    BOT1 = np.zeros((nfiles, N * Fc / Fa))
    TOP2 = np.zeros((nfiles, N * Fc / Fa))
    BOT2 = np.zeros((nfiles, N * Fc / Fa))
    i=0
    for x in l:
        (p1, p2) = findMaxPhase(x, total_t, Fa, Fc, False)
        if np.size(p1) > 1 or np.size(p2) > 1:
            print('Phase potentially zero for file = %s' % x)
        sq1 = sig.square(p1 + base_phase)
        sq2 = sig.square(p2 + base_phase)
        print('Phase 1 = %f' % p1)
        print('Phase 2 = %f' % p2)
        d = np.array(h5py.File(x)['PMT_DATA_8BIT'])
        s1 = sum(d[:,0:15])
        s2 = sum(d[:,16:])
        [T1, B1] = getTopBot(s1, square1, 1/Fa, Fc/2)
        [T2, B2] = getTopBox(s2, square2, 1/Fa, Fc/2)
        TOP1[i,:] = T1
        BOT1[i,:] = B1
        TOP2[i,:] = T2
        BOT2[i,:] = B2
    return (TOP1, BOT1, TOP2, BOT2)


def findMaxPhase(filename, Ta, Fs, Fc, plots, duty=0.5):
    '''
    FUnction which calculates and plots the maximum phase for the square wave.
This assumes that you are using the synchronization box, and generates a
square wave programmatically rather than using inverse transforms.

It will display plots of the phase for the user to verify that the phase
is correct and the same for both lasers.

[phase1, phase2] = findMaxPhase(filename, Ta, Fs, Fc, plots)
INPUTS:
   filename
   Ta          - Acq time in seconds
   Fs          - Sampling frequency
   Fc          - Chop frequency
   plots       - Boolean to display plots
    '''
    d = np.array(h5py.File(filename, 'r')['PMT_DATA_8BIT'])
    s1 = np.sum(d[:,6:9], 1)
    s2 = np.sum(d[:,22:25], 1)
    phases = np.linspace(0, 2 * np.pi * Fc * Ta, Ta * Fs)
    p = np.linspace(0, 2 * np.pi, 100)
#    PMT1 = np.zeros((100,1))
#    PMT2 = np.zeros((100,1))
    PMT1 = [np.sum((sig.square(phases + x, duty) + 1) * s1) / 2 for x in p]
    PMT2 = [np.sum((sig.square(phases + x, duty) + 1) * s2) / 2 for x in p]
    im1 = np.where(PMT1 == np.max(PMT1))
    phase1 = p[im1]
    im2 = np.where(PMT2 == np.max(PMT2))
    phase2 = p[im2]
    return (phase1, phase2)

def GetTopBot(sums, squares, dt, nyq):
    '''
%[TOP, BOT] = getTopBot(sums, squares, dt, nyq)
Function which takes the signal, aligned square wave, and time
differential and generates the top and bottom arrays.

(TOP, BOT) = GetTopBox(sums, squares, dt, nyq)

INPUTS:
sums       = Array of pmt readings.
squares    = Square wave.
dt         = time difference between readings. 1 / f
nyq        = Desired nyquist frequency (usually 1/2 of chop)
OUTPUTS:
TOP        = Points corresponding to the downsampled laser ON
BOT        = Poitns corresponding to the downsampled laser OFF
'''
#Old
#    su = sums * (squares + 1) / 2
#    sb = sums * (-squares + 1) / 2
#    [f, gu] = spec.spec(su, dt)
#    [f, gb] = spec.spec(sb, dt)
#    cut = np.where(abs(f) > nyq)
#    gu[cut] = 0
#    gb[cut] = 0
#    [t, su_c] = ispec(gu, f)
#    [t, sb_c] = ispec(gb, f)
#    #Do I need to put in a shift? Take a look at some data.
#    TOP = sig.resample(su_c, len(su_c) * nyq * 2 * dt)
#    BOT = sig.resample(sb_c, len(sb_c) * nyq * 2 * dt)
#New
    #Downsample the array
    [f, g] = spec.spec(sums)
    cut = np.where(abs(f) > nyq)
    g[cut] = 0
    [t, sums_ds] = spec.ispec(g, f)
    #Determine the duty cycle of the square wave.
    # (# of differences per sq. wave cycle) * (scaling factor 
    # for all cycles in
    # square array) * normalization + (default duty cycle)
    d = (np.sum(squares) * (10 / len(squares)) / 20.0) + 0.5
    d = round(d, 1)
    #Determine phase, in multiple of 2pi/10
    phs = round(10 * d - np.where(squares[0:10] == -1)[0][0])
    #Apply the square wave to the data.
    su = sums_ds * (squares + 1) / 2
    sd = sums_ds * (squares + 1) / 2
    if d == 0.5:
        TOP = sig.resample(su, len(su) * nyq * 2 * dt)
        BOT = sig.resample(sd, len(sd) * nyq * 2 * dt)
    else:
        on = d * 10
        off = (1 - d) * 10
        #Get rid of zerod elements
        su = su[(su != 0)]
        sd = sd[(sd != 0)]
        su = su.reshape(on, len(su) / on)
        sd = sd.reshape(off, len(sd)/off)
        #Weighted average along each axis.


