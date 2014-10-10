#This is a basic script to put together a bunch of response functions for a given directory. 
import numpy as np
import scipy
from scipy import io, signal, fftpack as fp
import glob
import os
import linear_response
import Gnuplot

l = glob.glob('*cm.mat')

#Get the linear response functions
#arrays = [(x['xcsmean_lifscaled'], x['ac_ch1_mn'], x['ac_ch2_mn']) for x in [io.loadmat(y) for y in l]]

#Load the arrays out of the files, grab xcsmean_lifscaled, ac_ch1_mn, and ac_ch2_mn, and get the linear 
#response functions. This one's a bit of a handful.
r_t = [linear_response.response(xx, yy, zz, 0.01, 0.04, 1e-6) for xx, yy, zz in [(x['xcsmean_lifscaled'][0,:], x['ac_ch1_mn'][0,:], x['ac_ch2_mn'][0,:]) for x in [io.loadmat(y) for y in glob.glob('*cm.mat')]]]

Hxy = []
Hyx = []
t = []
ct = 0

for x in r_t:
    Hxy.append(x[0])
    Hyx.append(x[1])
    t.append(x[2])

#Let's make some plots!
g1 = Gnuplot.Gnuplot(debug=1)
g2 = Gnuplot.Gnuplot(debug=1)

for x in r_t:
    Hxy = x[0]
    Hyx = x[1]
    t = x[2]
    fn = l[ct]
    num = '%0*d' % (3, ct)
    g1('set output "plots/Hxy_' + str(num) + '.png"')
    g2('set output "plots/Hyx_' + str(num) + '.png"')
    ct += 1
