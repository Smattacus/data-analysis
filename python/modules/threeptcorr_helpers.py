import ctypes
mylib = ctypes.cdll.LoadLibrary('/Users/smattingly/Programs/data-analysis/c/modules/p_xcorr/libxcorr.so')

import numpy as np

def getAvgCorr(data, delta_t):
    tau1 = c_int(delta_t)
    tau2 = c_int(delta_t)
    arrlarge = c_double * ((2 * delta_t - 1)**2)
    runavgarr = arrlarge()
    bigarr = arrlarge()
    for d1, d2, a1, a2 in zip(d['D1_scale'], d['D2_scale'], d['ave1'].transpose(), d['ave2'].transpose()):
        f1 = d1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        f2 = d2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        if a1 >= a2:
#            mylib.corr3_parallel(f1, f2, a1, a2, tau1, tau2, 1, bigarr)
            print(a1)
        else:
#            mylib.corr3_parallel(f2, f1, a2, a1, tau1, tau2, 1, bigarr)
            print(a2)
        runavgarr += bigarr
        d1
        a1
    runavgarr /= np.size(d['ave1'])
    return runavgarr
