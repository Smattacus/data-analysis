#Script for debugging the demodulation process.

import glob
import lifxcmeans
import flucsim
import numpy as np
import h5py
import scipy.signal as sig
from matplotlib.pyplot import *

dir = '/home/sean/Data/07-02-2016'

wnl = glob.glob(dir + '/*WhiteNoise_500kC*')
wnl.sort()

TB = lifxcmeans.genTPBoxFindPhaseOld(wnl[0], 3, 1e6, 1e5, duty=0.5,  plots=True)

(p1, p2) = lifxcmeans.findMaxPhaseViaSum(wnl[0], 3, 1e6, 1e5, plots=True)

(p1_old, p2_old) = lifxcmeans.findMaxPhase(wnl[0], 3, 1e6, 1e5, plots=True)

d = np.array(h5py.File(wnl[0], 'r')['PMT_DATA_8BIT'])

phasegen = np.linspace(0, 2 * np.pi * 1e5 * 3, 3 * 1e6, endpoint=False)

sq1 = sig.square(phasegen+ p1, 0.5)
sq2 = sig.square(phasegen + p2, 0.5)

sq1_old = sig.square(phasegen + p1_old, 0.5)
sq2_old = sig.square(phasegen + p2_old, 0.5)

s1 = np.sum(d[:,0:16], 1)
s2 = np.sum(d[:,16:], 1)


