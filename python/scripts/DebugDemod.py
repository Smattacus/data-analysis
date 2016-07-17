#Script for debugging the demodulation process.

import glob
import lifxcmeans
import flucsim
import numpy as np
import h5py
import scipy.signal as sig
from matplotlib.pyplot import *
import spec
import scikits.samplerate

dir = '/home/sean/data/SeanRunOne/Raws/04-03-15/Diode_P0Step/Dye_611p659794/'

fn = 'Fri_Apr_03_20_08_31_num5_2p5e-4torr_RP29WF_Diode_668p6114024_Dye611p659794.h5'
# Afaict, my phase max finding routine is working properly. So use it.

d = np.array(h5py.File(dir + fn, 'r')['PMT_DATA_8BIT'])

(p1, p2) = lifxcmeans.findMaxPhase(dir + fn, 8, 1e6, 1e5, plots=True, duty=0.5)

phases = lifxcmeans.genBasePhase(8, 1e6, 1e5)

sq1 = sig.square(p1 + phases)
sq2 = sig.square(p2 + phases)

s1 = np.sum(d[:, 0:16], 1)
s2 = np.sum(d[:, 16:], 1)

su1 = s1 * (sq1 + 1) / 2
sb1 = s1 * (-sq1 + 1) / 2

[f, gu1] = spec.spec(su1, 1e-6)
[f, gb1] = spec.spec(sb1, 1e-6)

cut = np.where(np.abs(f) > 5e4)

gu1[cut] = 0
gb1[cut] = 0

[t, su1_c] = spec.ispec(gu1, f[1] - f[0])
[t, sb1_c] = spec.ispec(gb1, f[1]- f[0])

top1 = resample(su1_c, 0.1, 'sinc_best')
bot1 = resampel(sb1_c, 0.1, 'sinc_best')