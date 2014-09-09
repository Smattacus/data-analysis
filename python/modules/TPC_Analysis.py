import numpy as np
import spec

def Get_TPC_Spec(tpc, dt, applywindow = True, window = 'gaussian',
        window_var=0.01, returnwin=False):
    """
    Calculates a windowed bispectrum of a three point correlation matrix.
    Returns the bispectrum, the frequency axis, and optionally the time and
    window data.
    INPUTS:
        tpc         
            - Three point correlation data. 2D numpy array.
        dt          
            - Time differential between elements in the TPC. float.
        applywindow=True 
            - Boolean. Whether to apply a gaussian window or not.
        window_var=0.01      
            - Float. Window size on gaussian window.
        returnwin=False
            - Boolean. Whether to return the window or not.

        OUTPUTS:
        f1
            - np array of frequency data along axis 1.
        f2
            - np array of frequency data along axis 2.
        bsp
            - 2D np array of (Windowed) bispectrum.
       win
            - 2D np array of window if returnwin=True

    """
#    This is an internal boolean. This will be added to the arguments should we
#   decide to use non-gaussian windows.
    window = 'gaussian'
    N1 = tpc.shape[0]
    N2 = tpc.shape[1]
    tc = window_var
    if applywindow:
        t1 = np.linspace(-(N1 - 1)/2, (N1 - 1)/2, N1) * dt
        t2 = np.linspace(-(N2 - 1)/2, (N2 - 1)/2, N2) * dt
        if window == 'gaussian':
            [Xg, Yg] = np.meshgrid(t1, t2)
            corr2win = np.exp(-(Xg / tc)**2 - (Yg / tc)**2 - (Xg * Yg) / (tc)**2)
        tpcwin = tpc * corr2win
        [f1, f2, g] = spec.spec2d(tpcwin, dt)
    if returnwin:
        return f1, f2, g, corr2win
    else:
        return f1, f2, g

def surface_tpc(mat, x, y, gridsize = (20, 20), imagename=''):
    """
        Quick helper routine to plot up a section of the bispectrum using
        Gnuplot. Will write to a file if imagename is given a value.
        mat
            - 2d np array. Array to be plotted.
        x
            -
"""
    return

