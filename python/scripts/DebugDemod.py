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

s1 = np.sum(d[:, 0:16], 1)
s2 = np.sum(d[:, 16:], 1)

# Now let's start looking at the demodulated data for one file, and
# figure out why it looks weird.

TB = np.array(TB)
T1 = TB[0, 0, :]
B1 = TB[0, 1, :]
T2 = TB[1, 0, :]
B2 = TB[1, 1, :]

# Try out gettopbot with the older code:
top1, bot1 = lifxcmeans.getTopBot(s1, sq1_old, 1e-6, 1e5/2, p1_old, 0.5)
top2, bot2 = lifxcmeans.getTopBot(s2, sq2_old, 1e-6, 1e5/2, p2_old, 0.5)

