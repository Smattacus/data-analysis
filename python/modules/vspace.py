#Some functions to streamline the process of 2D deconvolution. I've been
#writing too many scripts and I need moree on the fly flexibility.

#This is also going to be needed to streamline the process of averaging(?)
#along the time dimension of the 2d cross correlation / generating movies.
import numpy as np
import zeeman 
from scipy import interpolate as interp
import subprocess

def gen2dSpecs(x1, x2, x1dv, x2dv, retlines=False):
    '''
    Generates the two dimensional meshgrid of spectra to be used in
    deconvolution.
    
    INPUTS:
    x1      :   X axis wavelengths (Dye wavelengths)
    x2      :   Y axis wavelengths (Diode wavelengths)
    x1dv    :   Zeeman data value object for data array 1 (Zeeman.Dye_Values()
                or Zeeman.Diode_Values() output)
    x2dv    :   Zeeman data value object for data array 2.
    '''
    d1_glx, d1_wlpd, d1_lns = zeeman.ZeemanSpec_Padded(x1, np.array([1, 1]), (-1,),
            1e3, *x1dv)
    d2_glx, d2_wlpd, d2_lns = zeeman.ZeemanSpec_Padded(x2, np.array([1,1]), (-1,),
            1e3, *x2dv)
    [gx1, gx2] = np.meshgrid(d1_glx, d2_glx, indexing='xy')
    if retlines == True:
        return (gx1, gx2, d1_lns, d2_lns)
    else:
        return (gx1, gx2)

def getDataCoords(x1, x2):
    '''
    Generates the data coordinate array. This is assuming that x is a list of
    arrays, while y is just 21 values. This is used for the dye data which is
    spread over many wavelengths, while the diode is just spread of 21.
    '''
    data_coords = np.vstack([[(y, x[1]) for y in x[0]] for x in zip(x1, x2)])
    return data_coords

def genDataGrids(x1, x2, xcorrs, index, numfact = 1e-5, indexing='xy',
        method='linear', pad=0, fill_value=0):
    '''
        Function which generates an interpolated grid for a certain index of the 2
        dimensional cross correlation function. 
        INPUTS:
        x1      :   x dimension to meshgrid on. Assumed to be a list of arrays
                    - this describes the dye pattern of measurements.
        x2      :   y dimension to meshgrid on.
        xcorrs  :   Cross correlation matrix [N x N x (2 * t - 1)]
        index   :   which 2D xcorr plane to use (xc[:,:,index])
        numfact :   Factor to determine number of points to make a grid over.
                    Determined for each axis:
                    x1num   = (np.max(x1) - np.min(x1)) / 1e-6
                    x2num   = (np.max(x2) - np.min(x2)) / 1e-6
        indexing='xy'   :   meshgrid option for xy or ij indexing. Defaults to
                            xy.
        interp='linear' :   Interpolation for meshgrid. Defaults to 'linear'.
        OUTPUT:
        wlx     :   gridded x data points.
        wly     :   gridded y data points.
        gr      :   Gridded and interpolated data.
        data    :   data from xcorr[:,:,index]
        data_coords :   x and y coordinates corresponding to the xcorr array
                        elements.
    '''
    data = xcorrs[:,:, index].reshape(441,1) #Progresses along each x axis, then advances.
    dc = getDataCoords(x1, x2)
    x1d = np.linspace(np.min(dc[:,0]) - pad, np.max(dc[:,0]) + pad, num=(np.max(x1) -
        np.min(x1) + 2 * pad)/numfact)
    x2d = np.linspace(np.min(dc[:,1]) - pad, np.max(dc[:,1]) + pad, num=(np.max(x2) -
        np.min(x2) + 2 * pad)/numfact)
    wlx, wly = np.meshgrid(x1d, x2d, indexing=indexing)
    gr = interp.griddata(dc, data, (wlx, wly), method=method,
            fill_value=fill_value)
    return (wlx, wly, gr, data, dc, x1d, x2d)

def makeWebm(ds, fileprefix, sigfig):
    '''
        Calls the bash script make_webm.sh in order to automatically create a
        webm of PNGs that have been generated with fileprefix.

        The filename is assumed to have the form

        ${ds}${fileprefix}_%0${sigfig}.png

        INPUTS:
        ds          :   Directory string. Should have trailing slash.
        fileprexif  :   Prefix used to write the PNG files.
        sigfig      :   Numer of digits in the file numbers. I.e., 0045
                        corresponds to sigfig = 4. (aka written with %04d)
                        

    '''
    if ds[-1] != '/':
        ds += '/'
    subprocess.call(['/home/sean/programs/data-analysis/sh_scripts/make_webm.sh', ds + fileprefix, fileprefix + '.webm',
        str(sigfig)])
