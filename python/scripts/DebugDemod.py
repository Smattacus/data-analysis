#Script for debugging the demodulation process.

import glob
import lifxcmeans
import flucsim
import numpy as np
import h5py
import scipy.signal as sig
from matplotlib.pyplot import *

dir = '/home/sean/data/SeanRunOne/Raws/04-03-15/Diode_P0Step/Dye_611p659794/'

fn = 'Fri_Apr_03_20_08_31_num5_2p5e-4torr_RP29WF_Diode_668p6114024_Dye611p659794.h5'
# Afaict, my phase max finding routine is working properly. So use it.

d = np.array(h5py.File(dir + fn, 'r')['PMT_DATA_8BIT'])

(p1, p2) = lifxcmeans.findMaxPhase(dir + fn, 8, 1e6, 1e5, plots=True, duty=0.5)

