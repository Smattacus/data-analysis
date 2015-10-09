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

def genAllTPBoxFindPhase(path, Ta, Fs, Fc, duty=0.6, chan=0, ch1os=1, ch2os=1):
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
   path     : Path to a directory of h5 files.
   Ta       :   Total # of seconds of acquisition
   Fc       : Frequency of the square wave function. Generally, this is the
            laser chop frequency
   Fa       : Acquisition frequency.
   duty     : Duty cycle of square wave to use in signal matching. Empirically,
   .60 is the way to go.
   chan=0   :   Which PMT channels. to analyze. chan=0 (default) does both.
                    chan = 1 : PMT 1
                    chan = 2 : PMT 2
    ch1os   :   Displacement of channel 1 from the sumToTen() max. Check
                    figure for verification of correct value.
    ch2os   :   Displacement of channel 1 from the sumToTen() max. Check
                    figure for verification of correct value.
            NB: ch1os and ch2os default to 1, which means there is 1 sample of
                rising LIF signal before the max is hit.
OUTPUTS:
    Depending on chan argument, either (TOP1, BOT1), (TOP2, BOT2) or (TOP1,
    BOT1, TOP2, BOT2).
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
    TB = np.array([genTPBoxFindPhase(x, Ta, Fs, Fc, duty=duty, chan=chan,
            ch1os=ch1os,ch2os=ch2os) for x in l])
    if chan == 0:
        return (TB[:,0,0,:], TB[:,0,1,:], TB[:,1,0,:], TB[:,1,1,:])
    if chan == 1 or chan == 2:
        return (TB[:,0,:], TB[:,1,:])

def genTPBoxFindPhase(fn, Ta, Fs, Fc, plots=False, duty = 0.6, chan=0, ch1os =
        1, ch2os = 1):
    '''Generates the TOPs and BOTs of a single file for a given channel or both
    channels.

    INPUTS:
        fn      :   Filename of the relevant data.
        Ta      :   Total seconds of acquisition.
        Fs      :   Sampling speed of the PMT counting circuitry.
        Fc      :   Chop speed of the laser.
        plots   :   Whether to make plots of the phases. Defaults to false.
        duty=0.6:   Duty cycle of the square wave for alignment. Defaults to
                    0.6.
        chan=0  :   Which channel to get results for. Defaults to both channels
                    (chan = 0). chan = 1 corresponds to PMT1, while chan = 2
                    corresponds to PMT2.
        ch1os   :   Displacement of channel 1 from the sumToTen() max. Check
                    figure for verification of correct value.
        ch2os   :   Displacement of channel 1 from the sumToTen() max. Check
                    figure for verification of correct value.
                NB: ch1os and ch2os default to 1, which means there is 1 sample of
                rising LIF signal before the max is hit.
        OUTPUTS:
        if chan == 0:
            returns ([T1, B1], [T2, B2])
        if chan == 1:
            returns (T1, B1)
        if chan == 2:
            returns (T2, B2)
    '''
    (p1, p2) = findMaxPhaseViaSum(fn, Ta, Fs, Fc, plots, chan=chan, ch1os =
            ch1os, ch2os = ch2os) 
    d = np.array(h5py.File(fn, 'r')['PMT_DATA_8BIT'])
    phasegen = np.linspace(0, 2 * np.pi * Fc * Ta, Ta * Fs, endpoint = False)
    sq1 = sig.square(phasegen + p1, duty)
    sq2 = sig.square(phasegen + p2, duty)
    s1 = np.sum(d[:,0:16],1)
    s2 = np.sum(d[:,16:],1)
    if chan == 0:
        (T1, B1) = getTopBot(s1, sq1, 1 / Fs, Fc / 2, p1, duty)
        (T2, B2) = getTopBot(s2, sq2, 1/Fs, Fc / 2, p2, duty)
        return ([T1, B1], [T2, B2])
    if chan == 1:
        (T1, B1) = getTopBot(s1, sq1, 1/ Fs, Fc/2, p1, duty)
        return (T1, B1)
    if chan == 2:
        (T2, B2) = getTopBot(s2, sq2, 1/Fs, Fc/2, p2, duty)
        return (T2, B2)
    return 

def genTPBoxFindPhaseOld(fn, Ta, Fs, Fc, plots=False, duty = 0.6, chan=0):
    '''Generates the TOPs and BOTs of a single file for a given channel or both
    channels.

    INPUTS:
        fn      :   Filename of the relevant data.
        Ta      :   Total seconds of acquisition.
        Fs      :   Sampling speed of the PMT counting circuitry.
        Fc      :   Chop speed of the laser.
        plots   :   Whether to make plots of the phases. Defaults to false.
        duty=0.6:   Duty cycle of the square wave for alignment. Defaults to
                    0.6.
        chan=0  :   Which channel to get results for. Defaults to both channels
                    (chan = 0). chan = 1 corresponds to PMT1, while chan = 2
                    corresponds to PMT2.
                NB: ch1os and ch2os default to 1, which means there is 1 sample of
                rising LIF signal before the max is hit.
    '''
    (p1, p2) = findMaxPhase(fn, Ta, Fs, Fc, plots, duty=duty) 
    d = np.array(h5py.File(fn, 'r')['PMT_DATA_8BIT'])
    phasegen = np.linspace(0, 2 * np.pi * Fc * Ta, Ta * Fs, endpoint = False)
    sq1 = sig.square(phasegen + p1, duty)
    sq2 = sig.square(phasegen + p2, duty)
    s1 = np.sum(d[:,0:16],1)
    s2 = np.sum(d[:,16:],1)
    if chan == 0:
        (T1, B1) = getTopBot(s1, sq1, 1 / Fs, Fc / 2, p1, duty)
        (T2, B2) = geetTopBot(s2, sq2, 1/Fs, Fc / 2, p2, duty)
        return ([T1, B1], [T2, B2])
    if chan == 1:
        (T1, B1) = getTopBot(s1, sq1, 1/ Fs, Fc/2, p1, duty)
        return (T1, B1)
    if chan == 2:
        (T2, B2) = getTopBot(s2, sq2, 1/Fs, Fc/2, p2, duty)
        return (T2, B2)
    return 

def genTBChanPickPhase(filename, Ta, Fs, Fc, p1, p2, chans, plots=False, duty=0.5)
    ''''
    Generates a tuple of TOPS and BOTS corresponding 
    '''

def genTBChanSum(filename, Ta, Fs, Fc, plots=False, duty=0.6, s1range=(0,16),
        s2range=(16,32), ch1os=1, ch2os=1):
        '''
        Method for finding the TOPs and BOTs given a certain range of PMT data
        to run over. Instead of assuming that the two separate PMTs are from
        channels [0:16) and [16:32) 
        filename        :   Input file to read.
        Ta              :   Number of seconds data were taken over
        Fs              :   Sampling speed of data
        Fc              :   Chop frequency of laser
        plots=False     :   Generate plots? Defaults to False.
        duty=0.6        :   Duty cycle of square wave. Defaults to 60 % on.
        s1range=(0:16)  :   Indices for the first data set.
        s2range=(16:32) :   Indices for the second data set.

        If s2range is set to None, then the routine will find only one set
        of TOPS / BOTTOMS.
        
        The phase value used for s1range and s2range will be picked depending
        on which PMT the ranges fall in.
        '''
        d = np.array(h5py.File(filename, 'r')['PMT_DATA_8BIT'])
        s1 = np.sum(d[:,s1range[0]:s1range[1]],1)
        #Get the phase for each PMT.
        (p1, p2) = findMaxPhaseViaSum(filename, Ta, Fs, Fc, plots=False,
            plotname='temp.png', chan=0, ch1os=ch1os, ch2os=ch2os)
        phasegen = np.linspace(0, 2 * np.pi * Fc * Ta, Ta * Fs, endpoint =
            False)
        #If s2range is on PMT1, use phase 1
        if s2range[0] < 16 and s2range[1] << 16:
            p2 = p1
        #If s1range is on PMT2, use phase 2.
        if s1range[0] >= 16 and s1range[1] >= 16:
            p1 = p2
        if s2range == None:
             (T, B) = getTopBot(s1, sq, 1 / Fs, Fc / 2, p1, duty=0.6,
                 osample=10) 
             return (T, B)
        sq = sig.square(p1 + phasegen, duty) 
        s2 = np.sum(d[:,s2range[0]:s2range[1]],1)
        (T1, B1) = getTopBot(s1, sq, 1/Fs, Fc / 2, p1, duty = 0.6, osample=10)
        (T2, B2) = getTopBot(s2, sq, 1/Fs, Fc/2, p2, duty=0.6, osample=10)
        return (T1, B1, T2, B2)

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
    #No endpoint; else the spacing is all wonky.
    #phases = np.linspace(0, 2 * np.pi * Fc * Ta, Ta * Fs, endpoint=False)
    phases = genBasePhase(Ta, Fs, Fc)
    #101 points makes the ends both include zero.
    p = np.linspace(0, 2 * np.pi, 101, endpoint=True)
    PMT1 = [np.sum((sig.square(phases + x, duty) + 1) * s1) / 2 for x in p]
    PMT2 = [np.sum((sig.square(phases + x, duty) + 1) * s2) / 2 for x in p]
    if plots == True:
        pyplot.figure(1)
        pyplot.clf()
        pyplot.plot(p, PMT1, '-*')
        pyplot.xlabel('Phase offset')
        pyplot.ylabel('Sq wave xc')
        pyplot.title('PMT1 Square wave correlation')
        pyplot.figure(2)
        pyplot.clf()
        pyplot.plot(p, PMT2, '-*')
        pyplot.xlabel('Phase offset')
        pyplot.ylabel('Sq wave xc')
        pyplot.title('PMT2 Square Wave Correlation')
    im1 = np.where(PMT1 == np.max(PMT1))
    phase1 = p[im1]
    im2 = np.where(PMT2 == np.max(PMT2))
    phase2 = p[im2]
    if len(phase1) > 1:
        phase1 = phase1[0]
    if len(phase2) > 1:
        phase2 = phase2[0]
    return (phase1, phase2)

def findMaxPhaseViaSum(filename, Ta, Fs, Fc, plots=False,
plotname='temp.png', chan=0, ch1os=1, ch2os=1):
    '''
    Makes a histogram of the file using sumToTen(). A major caveat in this
    routine is that it ASSUMES that the data takes one sample to rise to the
    maximum LIF value. Most data taken at 10:1 Fs:Fc (1 MHz sample, 100 KHz
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
        ch1os=1     :   Offset of channel 1 relative to the max of PMT1.
                        Defaults to 1, which means there is 1 sample of rising
                        LIF before the max.
        ch2os=1     :   Offset of channel 2 relative to the max of PMT2.
                        Defaults to 1, which means is 1 sample of rising LIF
                        before the max.
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
    p1 = 2 * np.pi - (s1m-ch1os) * 2 * np.pi / (Fs / Fc)
    p2 = 2 * np.pi - (s2m-ch2os) * 2 * np.pi / (Fs / Fc)
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
    cut = np.where(abs(f) > nyq)
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
    on = duty * osample
    off = (1 - duty) * osample
    #Get rid of zerod elements and average over remaining columns
    TOP = np.mean(su.reshape(len(su) / osample, osample)[:, 0:on], 1)
    BOT = np.mean(sd.reshape(len(sd) / osample, osample)[:, on:],1)
    return (TOP, BOT)

def sumDataToTen(d1, d2, Ta, Fs, Fc):
    '''
    Routine for histogram PMT data which has already been summed into
    individual arrays.
    INPUTS:
        d1  :   
    '''
    s1rs = np.mean(d1.reshape(Ta * Fc, Fs / Fc), 0)
    s2rs = np.mean(d2.reshape(Ta * Fc, Fs / Fc), 0)
    return (s1rs, s2rs)

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
    s1rs = np.mean(s1.reshape(Ta * Fc, Fs / Fc), 0)
    s2rs = np.mean(s2.reshape(Ta * Fc, Fs / Fc), 0)
    return (s1rs, s2rs)

def plotHists(flist, Ta, Fs, Fc, target = '', chan=0):
    '''
        Generates histograms of the data files of all files in a given list.
        Writes PNG files to the target directory.

        The filename defaults to the input filename with extension changed to
        .png

        INPUTS:
        flist    :   List of files to output.
        Ta      :   Total acquisition time of each file.
        Fs      :   Sampling speed of each file.
        Fc      :   Chopping speed of the laser in the file.
        target  :   Target directory. Defaults to current directory.
        chan=0  :   Which channels to plot. Defaults to both.
                Options:
                chan = 0    : Plot both channels
                chan = 1    : Plot only PMT 1
                chan = 2    : Plot only PMT 2
    '''
    pyplot.figure(1)
    pyplot.clf()
    for x in flist:
        (s1rs, s2rs) = sumToTen(x, Ta, Fs, Fc)
        pyplot.clf()
        fo = target +  x.split('/')[-1].replace('.h5', '.png')
        if chan == 0 or chan == 1:
            pyplot.plot(s1rs)
        if chan == 0 or chan == 2:
            pyplot.plot(s2rs)
        pyplot.legend(['PMT1', 'PMT2'])
        pyplot.title('Plot for ' + fo)
        pyplot.savefig(fo)

def xcorru(x, y):
    '''
    Takes an unbiased cross correlation.
    '''
    N = len(x)
    xc = np.correlate(x, y, mode='full')
    xc = xc / np.append(np.linspace(1, N, N), np.linspace(N-1, 1, N-1))
    return xc

def genBasePhase(Ta, Fs, Fc):
    '''
    Routine which generates a base phase array by stacking individual arrays 
    generated with linspace that span some number of points from 0 to 2pi.
    '''
    Ntot = Ta * Fs
    Nperiods = int(Fc * Ta)
    single = np.linspace(0, 2 * np.pi, Fs / Fc, endpoint=False)
    return np.hstack([single.T for i in range(Nperiods)])

    
