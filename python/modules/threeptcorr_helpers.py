from ctypes import *
import ctypes

mylib = cdll.LoadLibrary('/Users/smattingly/Programs/data-analysis/c/modules/p_xcorr/libxcorr.so')

import numpy as np

def getAvgCorr(d, delta_t):
    #Runs the three point correlation using my C library. 
    #
    #INPUTS: 
    #   d    - Matlab dict containing "D1_scale", "D2_scale", "ave1", "ave2"
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
    avect = 0
    for x, y in zip(d['ave1'].transpose(), d['ave2'].transpose()):
        if x > y:
            avect += 1
    #Do a little magic with zip()
    for d1, d2, a1, a2 in zip(d['D1_scale'], d['D2_scale'],
            d['ave1'].transpose(), d['ave2'].transpose()):
        f1 = d1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        f2 = d2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        #Use the stronger signal as f1 in f1 * f1(t - tau1) * f2(t - tau2)
        if avect >= np.size(d['ave1'])/2:
            mylib.corr3_parallel(f1, f2, np.size(d1), np.size(d2), tau1, tau2, 1, bigarr)
        else:
            mylib.corr3_parallel(f2, f1, np.size(d2), np.size(d1), tau1, tau2, 1, bigarr)
        runavgarr +=  runavgarr + np.ctypeslib.as_array(bigarr)
    runavgarr /= np.size(d['ave1'])
    runavgarr = runavgarr.reshape(2 * delta_t - 1, 2 * delta_t - 1).transpose()
    return runavgarr
