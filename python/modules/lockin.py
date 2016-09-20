#File for reading out lock in data taken with Jorge's progrma.

import numpy as np
from matplotlib.pyplot import *
ion()
show()
import zeeman as z
import clean

def read_lockin(filename):
    d = np.loadtxt(filename, skiprows=1)
    wl = d[:,0]
    x = d[:,1]
    y = d[:,2]
    p = d[:,3]
    return (wl, x, y, p)

def background_subtract(wl, x, y, wl_thresh, p=False):
    i = np.where(wl < wl_thresh)
    x_bg = np.mean(x[i])
    y_bg = np.mean(y[i])
    x_bgs = x - x_bg
    y_bgs = y - y_bg
    if p==True:
        clf()
        plot(wl, x)
        plot(wl, y)
        plot(wl[i], x[i])
        plot(wl[i], y[i])
        xlabel('Wavelength')
        title('Real and Imag with BG estimate sections')
        legend(['Real', 'Imag'])
    return (x_bgs,  y_bgs)

def get_velocity(wl, x_bgs, y_bgs, wl_lthresh, wl_hthresh, p=False):
    #Set wl0 as the point where the real component crosses the x axis.
    d = np.abs(x_bgs)
    i = np.where((wl > wl_lthresh) ^ (wl > wl_hthresh))
    ic = np.where(min(d[i]) == d[i])
    if p == True:
        plot(wl, x_bgs)
        plot(wl, y_bgs)
        plot(wl[i], x_bgs[i])
        plot(wl[i], y_bgs[i])
        plot(wl[i][ic], x_bgs[i][ic], '*')
        plot(wl[i][ic], y_bgs[i][ic], '*')
        title('Velocity Selection area')
        legend(['real', 'imag', 'real mid', 'imag mid'])
    wl0 = wl[i][ic]
    v = z.doppv(wl, wl0)
    return v, wl0


def make_plots(filename, wl_bgthresh, wl_vthreshes):
    '''
    This is a helper function for making f0 and f1 plots.
    :param filename: Name of the file to read.
    :param wl_bgthresh: Upper limit for background estimation. Starts from left side.
    :param wl_vthreshes: tuple of (wavelength_left, wavelength_lelft) for 0 point estimation.
    :return:
    '''
    (wl, x, y, p) = read_lockin(filename)
    figure(1)
    clf()
    plot(wl, x)
    plot(wl, y)
    xlabel('Wavelength(nm)')
    title('Raw Data')
    legend(['Real',  'Imag'])
    figure(2)
    clf()
    (x_bgs, y_bgs) = background_subtract(wl, x, y, wl_bgthresh, p=True)
    figure(3)
    clf()
    zr = rotate_zero(x_bgs, y_bgs)
    plot(wl, np.real(zr))
    plot(wl, np.imag(zr))
    x_zr = np.real(zr)
    y_zr = np.imag(zr)
    title('Rotated ')
    figure(4)
    clf()
    (v, wl0) = get_velocity(wl, x_zr, y_zr, wl_vthreshes[0], wl_vthreshes[1], p=True)
    figure(5)
    clf()
    plot(v, x_zr)
    plot(v, y_zr)
    xlabel('Velocity (cm/s)')
    title('Real and Imag BG Subtracted vs Velocity.')
    legend(['Real', 'Imag'])
    return (wl, x, y, x_bgs, y_bgs, x_zr, y_zr, v, wl0)


def make_sweep_plots(filename, freqs, fignum=1, npoints = 2000, sweep=False):
    '''
    Makes plots for a sweep over antenna frequency, where p is the voltage for
    current sweep frequency.
    :param filename: File to read out.
    :param freqs: [low frequency, high frequency]
    :return: Tuple containing (wl, x, y, x_clean, y_clean, pha_clean)
    '''
    (wl, x, y, p) = read_lockin(filename)
    figure(fignum, figsize=(8,8))
    clf()
    subplot(221)
    f = gen_freq_axis(p, npoints, freqs[0], freqs[1])
    plot(f, x,'.')
    plot(f, y, '.')
    plot(f, np.abs(x + 1j * y), '.')
    xlabel('Sweep voltage (V)')
    ylabel('Lock in signal')
    title('Raw Sweep Data')
    legend(['Real', 'Imag', 'Rvec'])
    subplot(222)
    if sweep == False:
        plot(f, np.angle(x + 1j * y))
        xlabel('Sweep Voltage (V)')
        ylabel('Phase (rad)')
        title('Raw Phase Data')
    elif sweep == True:
        plot(p)
        plot(f * np.max(p) / freqs[1])
        xlabel('Sweep #')
        ylabel('Sweep Ramp Voltage')
        legend(['Raw Voltage', 'Constructed Sweep'])
        title('Constructed Sweep and Actual Sweep')
    subplot(223)
    (fc, xc) = clean.clean(f, y)
    (fc, yc) = clean.clean(f, x)
    plot(fc, xc)
    plot(fc, yc)
    xlabel('Frequency')
    ylabel('Lockin response(v)')
    title('Cleaned Lock in Response')
    legend(['Real', 'Imag'])
    subplot(224)
    phac = np.angle(xc + 1j * yc)
    plot(fc, np.unwrap(phac))
    xlabel('Frequency')
    ylabel('Phase (rad)')
    title('Cleaned Phase')
    return(x, y, p, xc, yc, fc)




def rotate_zero(x, y):
    '''
    Rotate according to x + 1j * y so that sum (x + 1j * y) = 0.
    :param x: Real component
    :param y: Imaginary component
    :return:
    '''
    z = x + 1j * y
    zs = np.sum(z)
    zsa = np.angle(zs)
    zr = z / np.exp(1j * zsa)
    return zr

def gen_freq_axis(sw_v, npoints, fmin, fmax):
    '''
    Generates a frequency axis based on the sweep voltage. This function generates
    a linear range of frequencies to remove the noise which has been introduced by
    measuring the voltage.

    :param sw_v: sweep voltage array.
    :param npoints: Number of points per sweep, after which it repeats.
    :param fmin: Minimum sweep frequency.
    :param fmax: Maximum sweep frequency.
    :return:
    '''
    #Find the point where the sweep drops back down to fmin.
    sd = np.diff(sw_v)
    j = np.where(sd == np.min(sd))[0]
    if j > npoints:
        j = j % npoints
    #We want to create an array such that it starts at the same point in
    #the sweep as the sw_v array
    flarge = np.tile(np.linspace(fmin, fmax, npoints), np.ceil(sw_v.shape[0] / npoints) + 1)
    f = flarge[npoints - (j+1):npoints - (j+1) + sw_v.shape[0]]
    return f
