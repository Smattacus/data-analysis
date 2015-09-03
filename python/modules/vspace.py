#Some functions to streamline the process of 2D deconvolution. I've been
#writing too many scripts and I need moree on the fly flexibility.

#This is also going to be needed to streamline the process of averaging(?)
#along the time dimension of the 2d cross correlation / generating movies.
import numpy as np
import Zeeman as Z

def gen2dSpecs(x1, x2, x1dv, x2dv):
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
    d1_glx, d1_wlpd, d1_lns = Z.ZeemanSpec_Padded(x1, np.array([1, 1]), (-1,),
            1e3, *x1d)
    d2_glx, d2_wlpd, d2_lns = Z.ZeemanSpec_Padded(x2, np.array([1,1]), (-1,),
            1e3, *x2d)
    [gx1, gx2] = np.meshgrid(d1_glx, d2_glx, indexing='xy')
    return(gx1, gx2)

def genDataGrids(x1, x2, xcorrs, index, indexing='xy', interp='linear'):
    '''
        Function which generates a grid for a certain index of the 2
        dimensional cross correlation function. 
        INPUTS:
        x1      :   x dimension to meshgrid on.
        x2      :   y dimension to meshgrid on.
        xcorrs  :   Cross correlation matrix [N x N x 2 * t - 1]
        index   :   which 2D xcorr plane to use (xc[:,:,index])
        indexing='xy'   :   meshgrid option for xy or ij indexing. Defaults to
                            xy.
        interp='linear' :   Interpolation for meshgrid. Defaults to 'linear'.
    '''
    data = xcorrrs[:,:, index].reshape(441,1) #Progresses along each x axis, then advances.
    gr = 
    

