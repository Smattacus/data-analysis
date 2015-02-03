#
#his is a file containing major functions from MATLAB used for the purpose of
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
import matplotlib as mpl
from matplotlib import pyplot

def genAllTPBoxFindPhase(path, Fc, Fa, duty):
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
   path     - Path to a directory of h5 files.
   Fc       - Frequency of the square wave function. Generally, this is the
            laser chop frequency
   Fa       - Acquisition frequency.
   duty     - Duty cycle of square wave to use in signal matching. Empirically,
   .60 is the way to go.

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
    base_phase = np.linspace(0 , 2 * np.pi * total_t * Fc, N, endpoint=False)
    TOP1 = np.zeros((nfiles, N * Fc / Fa))
    BOT1 = np.zeros((nfiles, N * Fc / Fa))
    TOP2 = np.zeros((nfiles, N * Fc / Fa))
    BOT2 = np.zeros((nfiles, N * Fc / Fa))
    i=0
    for x in l:
        (p1, p2) = findMaxPhaseViaSum(x, total_t, Fa, Fc, False)
        if np.size(p1) > 1 or np.size(p2) > 1:
            print('Phase potentially zero for file = %s' % x)
        sq1 = sig.square(p1 + base_phase)
        sq2 = sig.square(p2 + base_phase)
        print('Phase 1 = %f' % p1)
        print('Phase 2 = %f' % p2)
        d = np.array(h5py.File(x)['PMT_DATA_8BIT'])
        s1 = np.sum(d[:,0:15], 1)
        s2 = np.sum(d[:,16:], 1)
        [T1, B1] = getTopBot(s1, sq1, 1/Fa, Fc/2, p1, duty, 10)
        [T2, B2] = getTopBot(s2, sq2, 1/Fa, Fc/2, p2, duty, 10)
        TOP1[i,:] = T1
        BOT1[i,:] = B1
        TOP2[i,:] = T2
        BOT2[i,:] = B2
    return (TOP1, BOT1, TOP2, BOT2)


def findMaxPhase(filename, Ta, Fs, Fc, plots, duty=0.5):
    '''
    Function which calculates and plots the maximum phase for the square wave.
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
    s1 = np.sum(d[:,0:16], 1)
    s2 = np.sum(d[:,16:], 1)
    phases = np.linspace(0, 2 * np.pi * Fc * Ta, Ta * Fs, endpoint=False)
    #No endpoint; else the spacing is all wonky.
    p = np.linspace(2 * np.pi / 100, 2 * np.pi, 100, endpoint=True)
    PMT1 = [np.sum((sig.square(phases + x, duty) + 1) * s1) / 2 for x in p]
    PMT2 = [np.sum((sig.square(phases + x, duty) + 1) * s2) / 2 for x in p]
    im1 = np.where(PMT1 == np.max(PMT1))
    phase1 = p[im1]
    im2 = np.where(PMT2 == np.max(PMT2))
    phase2 = p[im2]
    return (phase1, phase2)

def findMaxPhaseViaSum(filename, Ta, Fs, Fc, plots=False,
plotname='temp.png', chan=0, ch1os=1, ch2os=1):
    '''
    Makes a histogram of the file using sumToTen(). A major caveat in this
    routine is that it ASSUMES that the data takes one sample to rise to the
    maximum LIF value - most data taken at 10:1 Fs:Fc (1 MHz sample, 100 KHz
    chop) seems to follow this.

    This routine returns a phase which is a multiple of 2 * pi / 10 based on
    these assumptions.

    Moreover, the option exists to output .png plots with which one may verify
    that the data follows these assumptions. It is recommended to do so for
    large data runs.

    INPUTS:
        filename    :   Input file to find the phase of.
        Ta          :   Total sampling time in seconds.
        Fs          :   Sampling speed in Hz
        Fc          :   Laser chopping speed 100 KHz
        plots=False :   Whether to make verification plots to file.
        plotname='temp.png'     : Filename for the verification plots.
        chan=0      : Channels to generate plots for. 0 : Both. 1: PMT1. 2:
                        PMT2.
    OUTPUTS:
        (p1,p2)          :   Phases for (PMT1, PMT2)
    '''
    (s1rss, s2rss) = sumToTen(filename, Ta, Fs, Fc)
    if plots==True:
        pyplot.figure(1)
        pyplot.clf()
        if chan == 0:
            pyplot.plot(s1rss)
            pyplot.plot(s2rss)
            pyplot.legend(['PMT1', 'PMT2'])
            pyplot.savefig(plotname)
        elif chan == 1:
            pyplot.plot(s1rss)
            pyplot.legend(['PMT1'])
            pyplot.savefig(plotname)
        elif chan == 2:
            pyplot.plot(s2rss)
            pyplot.legend(['PMT2'])
            pyplot.savefig(plotname)
        else:
            print('Invalid chan value! LIF_XCMeans.findMaxPhaseViaSum called'
            ' in error.')
    s1m = np.where(s1rss == np.max(s1rss))[0][0]
    s2m = np.where(s2rss == np.max(s2rss))[0][0]
    #Assume a one element offset
    p1 = 2 * np.pi - (s1m-1) * 2 * np.pi / 10
    p2 = 2 * np.pi - (s2m-1) * 2 * np.pi / 10
    return (p1, p2)

def getTopBot(sums, squares, dt, nyq, phase, duty, osample=10):
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
phase       = Phase of the input square wave.
duty        = Duty Cycle of the input square wave. Assumed to be .1 precision.
osample     = Oversampling rate. Assumed to be 10.
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
    [f, g] = spec.spec(sums, dt)
    cut = np.where(abs(f) > 2 * nyq)
    g[cut] = 0
    [t, sums_ds] = spec.ispec(g, f)
    sums_ds = np.abs(sums_ds)
    #Determine the duty cycle of the square wave.
    # (# of differences per sq. wave cycle) * (scaling factor 
    # for all cycles in
    # square array) * normalization + (default duty cycle)
    #Determine phase, in multiple of 2pi/10
    phs = int(round(phase / ( 2 * np.pi / osample))) 
    #Apply the square waves.
    su = sums_ds * (squares + 1) / 2
    sd = sums_ds * (-squares + 1) / 2
    #Roll the incomplete cycles from the end and beginning
    #together, then eliminate them.
    su = np.roll(su, phs)[osample:]
    sd = np.roll(sd, phs)[osample:]
    if duty == 0.5:
        TOP = sig.resample(su, len(su) * nyq * 2 * dt)
        BOT = sig.resample(sd, len(sd) * nyq * 2 * dt)
    else:
        on = duty * osample
        off = (1 - duty) * osample
        #Get rid of zerod elements and average over remaining columns
        TOP = np.mean(su.reshape(len(su) / osample, osample)[:, 0:on], 1)
        BOT = np.mean(sd.reshape(len(sd) / osample, osample)[:, on:],1)
    return (TOP, BOT)


def sumToTen(filename, Ta, Fs, Fc):
    '''
    Takes an input data file and returns 2 10 element arrays containing the
    data histogram in Fs / Fc bins.
    INPUTS:
        filename    :   Name of the file to generate histogram of.
        Ta          :   Total acquisition time in the file.
        Fs          :   Sampling speed in Hz
        Fc          :   Laser chopping speed in Hz
    OUTPUTS
        s1rs        :   Histogram for PMT1
        s2rs        :   Histogram for PMT2
    '''
    d = np.array(h5py.File(filename, 'r')['PMT_DATA_8BIT'])
    s1 = np.sum(d[:, 0:16], 1)
    s2 = np.sum(d[:,16:], 1)
    s1rs = np.sum(s1.reshape(Ta * Fc, Fs / Fc), 0)
    s2rs = np.sum(s2.reshape(Ta * Fc, Fs / Fc), 0)
    return (s1rs, s2rs)

def plotDirectorySums(flist, Ta, Fs, Fc, target = ''):
    '''
        Generates histograms of the data files of all files in a directory.
        Writes PNG files to the target directory.

        The filename defaults to the input filename with extension changed to
        .png

        INPUTS:
        flist    :   List of files to output.
        Ta      :   Total acquisition time of each file.
        Fs      :   Sampling speed of each file.
        Fc      :   Chopping speed of the laser in the file.
        target  :   Target directory. Defaults to current directory.
    '''
    figure(1)
    clf()
    for x in flist:
        (s1rs, s2rs) = sumToTen(f, ta, Fs, Fc)
        clf()
        fo = f.replace('.h5', '.png')
        plot(PMT1)
        plot(PMT2)
        legend(['PMT1', 'PMT2'])
        title('Plot for ' + fo)
        savefig(fo)


