#!/usr/bin/python2

#Import some packages. Mainly the scipy stuff.
#Consider being more selective with imports in the future, for speed of course.

import numpy as np
import matplotlib as mp
import os
import glob
import scipy as sp
from scipy import io

print("Hello world!")

N = 800000

os.chdir('/home/sean/data/5-28-2014/XCMeans/')

#Create a 2D array for the normalized cross correlations to try and recreate
#Fred's broadening peaks.
l = np.sort(glob.glob("*.mat"))

#Read all the mat files, grab the xcsmean_lifscaled array, take the middle 5000
#points. Isn't Python neat?
xcorr = [io.loadmat(x)["xcsmean_lifscaled"][0, ((N-1) - 2048):((N-1)+2049)] for x in l]
xcorr = np.array(xcorr)

#Note that the indices go from 0 - N+1 on slicing; they do 
#not include the final element. This is due to indexing from zero.

#Initialize the time array.
t =np.linspace(-(N-1), N-1, 2 * N - 1) / 1e5 
t = t[(N-1) - 2048: (N-1) + 2049]



