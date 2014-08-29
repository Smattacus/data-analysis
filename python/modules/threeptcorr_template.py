#!/opt/Python-2.7/bin/python2
#Path from sys.executable of python2 after loading module python27

#################################
#
#This is a basic Python script to run the parallel xcorr on all the .mat files
#in a given folder. A variable gives the script writer the flexibility to
#choose variables the way s/he wants. 
#
#If the variable is more than one string element, the script will assume that
#it is a list of filenames and will instead iterate through that.
#
#One can also pass input variables using typical argv calling conventions with
#script call:
#
#threeptcorr.py <searchstring> <path>
#
#################################

import sys
import numpy as np
import scipy as sp
from scipy import io
import os
import string
import glob

import threeptcorr_helpers as tpc 

####################################
#
#Input variables.
#
####################################

#If this is a list, the script will check it and read the files given in the
#list. Otherwise it will run a ls command and get the files from the fs.
if np.size(sys.argv) == 1:
    ################CHANGE HERE#####################
    search = '*DiffArrays.mat'
    p = '/Users/smattingly/Data/5-28-2014/XCMeans/' 
elif np.size(sys.argv) == 3:
    search = sys.argv[1]
    p = sys.argv[2]
    print("search read in is: " + search)
    print("Input path is: " + p)
else:
    print("There are " + str(np.size(sys.argv)) + " inputs.")
    for x in sys.argv:
        print(x)


#If it is a filename, glob.glob(search) will just return that filename.
#If there is only one element in the search list, search by that element name.
if np.size(search) == 1:
    os.chdir(p)
    fl = glob.glob(search)
else:
    os.chdir(p)
    fl = search


####################################

print("Starting corr loops.")
print("np.size = " + str(np.size(fl)))
print(fl)
print(os.getcwd())
#TAU VALUE
delta_t = 1000
for x in fl:
    data = io.loadmat(x)
    threepoint = tpc.getAvgCorr(data, delta_t)
    out = string.replace(x, '.mat', '_3corr.mat')
    print("Will save to " + out)
    io.savemat(out, {'tpcorr': threepoint})
