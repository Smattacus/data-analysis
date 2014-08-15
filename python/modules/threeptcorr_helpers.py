from ctypes import *
import ctypes

mylib = cdll.LoadLibrary('/Users/smattingly/Programs/data-analysis/c/modules/p_xcorr/libxcorr.so')

import numpy as np

def getAvgCorr(data, delta_t):
    #Runs the three point correlation using my C library. 
    #
    #INPUTS: 
    #   data    - Matlab dict containing "D1_scale", "D2_scale", "ave1", "ave2"
    #   delta_t - tau value to correlate by. Matrix returned will be
    #               [2 * delta_t - 1] x [2 * delta_t - 1]
    #
    #OUTPUTS    - 
    #   runavgarr - Pointwise averaged three point correlation function.
    #
    #
    
    tau1 = c_int(delta_t)
    tau2 = c_int(delta_t)
    arrlarge = c_double * ((2 * delta_t - 1)**2)
    runavgarr = np.zeros(( 2 * delta_t - 1) **2)
    bigarr = arrlarge()
    #Do a little magic with zip()
    for d1, d2, a1, a2 in zip(data['D1_scale'], data['D2_scale'],
            data['ave1'].transpose(), data['ave2'].transpose()):
        f1 = d1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        f2 = d2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        #Use the stronger signal as f1 in f1 * f1(t - tau1) * f2(t - tau2)
        if a1 >= a2:
            mylib.corr3_parallel(f1, f2, np.size(d1), np.size(d2), tau1, tau2, 1, bigarr)
        else:
            mylib.corr3_parallel(f2, f1, np.size(d2), np.size(d1), tau1, tau2, 1, bigarr)
        runavgarr +=  runavgarr + np.ctypeslib.as_array(bigarr)
    runavgarr /= np.size(data['ave1'])
    runavgarr = np.ctypeslib.as_array(runavgarr).reshape(tau1,
            tau2).transpose()
    return runavgarr
