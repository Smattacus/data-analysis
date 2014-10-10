#This is a basic script to put together a bunch of response functions for a given directory. 
import numpy as np
import scipy
from scipy import io, signal, fftpack as fp
import glob
import os
import linear_response
import Gnuplot
from numpy import real


#This is to be used as a script from any directory.
l = glob.glob('*cm.mat')
l.sort()

#Get the linear response functions
#arrays = [(x['xcsmean_lifscaled'], x['ac_ch1_mn'], x['ac_ch2_mn']) for x in [io.loadmat(y) for y in l]]

#Load the arrays out of the files, grab xcsmean_lifscaled, ac_ch1_mn, and ac_ch2_mn, and get the linear 
#response functions. This one's a bit of a handful.
r_t = [linear_response.response(xx, yy, zz, 0.01, 0.04, 1e-6) for xx, yy, zz in
        [(x['xcsmean_lifscaled'][0,:], x['ac_ch1_mn'][0,:],
            x['ac_ch2_mn'][0,:]) for x in [io.loadmat(y) for y in l]]]

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

g1('set terminal png')
g2('set terminal png')

xymax = 0
yxmax = 0
xymin = 0
yxmin = 0
for x in r_t:
    xymax = np.max(real(x[0])) if np.max(real(x[0])) >= xymax else xymax 
    yxmax = np.max(real(x[1])) if np.max(real(x[1])) >= yxmax else yxmax
    xymin = np.min(real(x[0])) if np.min(real(x[0])) <= xymin else xymin
    yxmin = np.min(real(x[1])) if np.min(real(x[1])) <= yxmin else yxmin

g1('set yrange [' + str(xymin) + ':' + str(xymax) + ']')
g2('set yrange [' + str(yxmin) + ':' + str(yxmax) + ']')


for x in r_t:
    Hxy = x[0]
    Hyx = x[1]
    t = x[2]
    tc = np.where(t < 0.0015)
    fn = l[ct]
    num = '%0*d' % (3, ct)
    fncm = fn.replace('.mat', '')
    g1('set output "plots/Hxy_' + str(num) + '.png"')
    g2('set output "plots/Hyx_' + str(num) + '.png"')
    ct += 1
    g1.title('Linear Response of X to Y at ' + str(fncm))
    g1.xlabel('Time (s)')
    g1.ylabel('Response')
    g2.title('Linear Response of Y to X at ' + str(fncm))
    g2.xlabel('Time (s)')
    g2.ylabel('Response')
    g1.plot(Gnuplot.Data(t[tc], np.real(Hxy[tc]), binary=0))
    g2.plot(Gnuplot.Data(t[tc], np.real(Hyx[tc]), binary=0))
