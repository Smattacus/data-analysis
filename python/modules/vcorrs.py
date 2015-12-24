#This is a file defining classes for handling sets of velocity correlation.
#Classes are for the purpose of making it easier to handle the corresponding
#sets of xcorrs, velocities, wavelengths, data paths, spectra, windowing
#functions, etc.
#

import numpy as np, scipy as sp, scipy.io as io, scipy.interpolate as interp
from matplotlib import mlab, pyplot as plt

import sroh, spec, clean, zeeman as z

cdye = z.getAmpCents(-1, 1e3, *(z.Dye_Values()[0:4]))[1][0]
cdi = z.getAmpCents(-1, 1e3, *(z.Diode_Values()[0:4]))[1][0]

class cls_vcorrset:
   

    def __init__(self, diode_step, tpc = False): 
        
        #Establish the list of .mat files for that day.
        ddict = sroh.getDataDict() 
        ds = 'Diode_P%dStep' % diode_step if diode_step >= 0 else 'Diode_N%dStep' % abs(diode_step)
        self.dfl = ddict[ds]
        self.diode_step = diode_step
        
        self.diode_wl = np.array(sroh.getDiodePoints()[diode_step + 10])

        self.dye_wl = np.array([float(x.split('/')[-1].split('_')[-1].split(
            '.mat')[0].replace('p','.')) for x in self.dfl])
        
        relindex = sroh.getDyeIVDFRels()[diode_step + 10]

        self.shift = sroh.getRels()[relindex]

        self.di_ivdf = sroh.getDiIVDFFiles()[relindex]
        self.dye_ivdf = sroh.getDyeIVDFFiles()[relindex]

        self.dye_v = z.doppv(self.dye_wl, cdye)
        self.diode_v = z.doppv(self.diode_wl, cdi)

        return
        
    def loadMats(self, tc = 1, tstep = 1e-5, noscale = False):
        '''
            Populates a 21 x (2 * N - 1) self.xcorrs numeric array by reading the .mat files
            given by self.dfl.

            This also initializes the self.trim_xcorrs, which is used by all
            the major data routines. self.xcorrs is used to hold the larger
            xcorr.

            Optional inputs:
            tc = 1      : (float) Gives the clipping of the cross correlation
            to keep in the array. In seconds. 
        '''
        arrs = []
        lif_av1 = []
        lif_av2 = []
        for x in self.dfl:
            print('Reading from ' + x)
            d = io.loadmat(x)
            if noscale == True:
                xc = d['xcmean'][0]
            if noscale == False:
                xc = d['xcsmean_lifscaled'][0]
            lif_av1.append(d['ave1'][0])
            lif_av2.append(d['ave2'][0])
            N = (xc.shape[0] + 1) / 2
            t = np.linspace(-(N-1), N-1, 2 * N - 1) * tstep
            arrs.append(xc[np.where(np.abs(t) < tc)])
        self.xcorrs = np.vstack(arrs)
        self.trim_xcorrs = self.xcorrs
        self.t = t[np.where(np.abs(t) < tc)]
        self.trim_t = self.t
        self.lif_avgs1 = lif_av1
        self.lif_avgs2 = lif_av2
        
    def loadNpy(self, prefix = '', tstep=1e-5):
        '''
            Populates the self.xcorrs variable by reading an npy file.
            INPUTS:
            filename    :   (string) .npy file to load using np.load().
        '''
        fn = ('obj_' + prefix + '_Diode_P%dStep.npz' % self.diode_step if self.diode_step >= 0 else
        'obj_' + prefix + '_Diode_N%dStep.npz' % abs(self.diode_step))
        read = np.load(fn)
        self.xcorrs = read['arr_0']
        self.lif_avgs1 = read['arr_1']
        self.lif_avgs2 = read['arr_2']
        self.N = (self.xcorrs.shape[1] + 1) / 2
        self.t = np.linspace(-(self.N - 1) , self.N - 1, 2 * self.N - 1) * tstep
        return

    def saveNpy(self, prefix = ''):
        '''
            Saves the xcorr data to a filename (for faster object construction
            after the initial run).

        '''
        fn = ('obj_' + prefix + '_Diode_P%dStep.npz' % self.diode_step if self.diode_step >= 0 else 
         'obj_' + prefix + '_Diode_N%dStep.npz' % abs(self.diode_step))
        np.savez(fn, self.xcorrs, self.lif_avgs1, self.lif_avgs2)
        return

    def trimXcorrs(self, tc = 0.1, tstep=1e-5):
        '''
            This function trims the object's xcorrs down to a certain size.

            INPUTS:
            tc      :       (float) Time (s) to cut down the self.t and
                            self.xcorrs arrays to.
            tstep   :       Time step.
        '''
        self.trim_xcorrs = self.xcorrs[:, np.where(np.abs(self.t) < tc)[0]]
        self.trim_t = self.t[np.where(np.abs(self.t) < tc)]

    def setCleanedIVDFS(self, fwin = 1500):
        '''This routine takes the self.di_ivdf and self.dye_ivdf files and
        generates a cleaned (fft - zerod) IVDF in self.di_ivdf_clean and
        self.dye_ivdf_clean.
        
        Optional Inputs:
        fwin = 1500     :   Wavenumber outside which to set the spec to 0.

        '''
        #First the dye.
        d = np.loadtxt(self.dye_ivdf, skiprows=1)
        R = np.sqrt(d[:,2]**2 + d[:,3]**2)
        wl = d[:,0]
        (wlc, rc) = clean.clean(wl, R)
        [f, g] = spec.spec(rc, wlc[1] - wlc[0])
        g[np.where(np.abs(f) > fwin)] = 0
        [t, rcf] = spec.ispec(g, f[1] - f[0])
        self.dye_ivdf_c = (wlc, rcf)
        #Then the diode.
        if self.di_ivdf == '':
            #The Diode file is missing.
            return
        else:
            d = np.loadtxt(self.dye_ivdf, skiprows=1)
            R = np.sqrt(d[:,2]**2 + d[:,3]**2)
            wl = d[:,0]
            (wlc, rc) = clean.clean(wl, R)
            [f, g] = spec.spec(rc, wlc[1] - wlc[0])
            g[np.where(np.abs(f) > fwin)] = 0
            [t, rcf] = spec.ispec(g, f[1] - f[0])
            self.diode_ivdf_c = (wlc, rcf)
        return
        
    def dualLIFPlot(self, y = None):
        '''
            Makes a vertical LIF photon / s rate plot using internal lif
            averages.
        '''
        if y is None:
            y = self.dye_v - self.diode_v
        l1 = plt.plot(np.array([np.mean(x) for x in self.lif_avgs1]), y,
                color='red', label='PMT1')
        plt.xticks(color='red')
        plt.xlabel('ave photons / demod sample (photon kHz)', color='red')
        plt.title('PMT1 & 2 LIF Signal', y=1.08)
        plt.ylim(np.min(y), np.max(y))
        plt.twiny()
        l2 = plt.plot(np.array([np.mean(x) for x in self.lif_avgs2]), y,
                color='red', label='PMT2')


    def xcSpec(self, tw = 0.01):
        '''
            Takes the spectrum of the xcorr array - along the long (time) axis.
            Uses the untrimmed array for the actual spectrum.

            INPUTS:
            tw     :       (float) Width of Gaussian window to use (if any).
        '''
        win = np.exp(-(self.t / tw)**2/2)
        #We have to make a grid of windowing functions.
        [xxwin, yy] = np.meshgrid(win, range(self.xcorrs.shape[0]))
        #xxwin is now a 2d array which has each row corresponding to a
        #windowing function, ready to be pointwise multiplied by xcorrs.
        [self.f, self.g] = spec.spec(xxwin * self.xcorrs, self.t[1] -
                self.t[0], axis=1)
        return

    def setGrids(self, y = None):
        '''
        Creates grid coordinates self.tt and self.vv based on the current t and
        v structures. 
        
        INPUTS:
        y = 0       : Alternate y axis to use. Defaults to self.dye_v - self.diode_v.
        '''
        if y is None:
            y = self.dye_v - self.diode_v
        self.tt, self.vv = np.meshgrid(self.trim_t, y)
        Ns = self.tt.shape[0] * self.tt.shape[1]
        self.zc = np.vstack([np.reshape(self.tt, Ns), np.reshape(self.vv, Ns),
            np.reshape(self.trim_xcorrs, Ns)])
        return

    def xcMesh(self, y = None,  tc = 0.01,
            ngridx = 100, ngridy = 100, **kwargs):
        '''
            Makes a surface mesh of the xcorr data.
            INPUTS:
            fig         :   (figure object) Axes to use.
            y = 0       :   (float) Y axis.Recommend self.dye_v - self.diode_v
            tc = 0.01   :   (float) Point to cut time axis.
            ngridx = 100    :   # of interpolation points on x axis.
            ngridy = 100    :   # of interpolation points on y axis.
        '''
        if y is None:
            y = self.dye_v - self.diode_v
        self.setGrids(y)
        xi = np.linspace(np.min(self.trim_t), np.max(self.trim_t), ngridx)
        yi = np.linspace(np.min(y), np.max(y), ngridy)
        xx, yy = np.meshgrid(xi, yi)
        #The main problem is that griddata wants equal length x and y arrays.
        #zg = mlab.griddata(self.zc[0], self.zc[1], self.zc[2], xi, yi, **kwargs)
        zg = interp.griddata((self.zc[0], self.zc[1]), self.zc[2], (xx, yy),
                **kwargs)
        plt.contour(xi, yi, zg, 15, linewidths=0.5, colors='k')
        plt.contourf(xi, yi, zg, 15)
        return (xi, yi, zg)

    def xcWinMesh(self, tl, th, y = None, ngridx = None, ngridy = None,
            **kwargs):
        '''
            Makes a surface mesh of a user selected range of xcorr data.

            INPUTS:
            tl      :       (float) Start time for the xcorr array.
            th      :       (float) End time for the xcorr array.
            y (None):       array of points to use for y axis. Defaults to
                            self.dye_v - self.diode_v
            ngridx   :      Number of grid points in x direction. Defaults to
                            10 * len(xcorr[tw])
            ngridy   :      Number of grid points in y diection. Defaults to
                            10 * len(y)
            **kwargs :      kwargs fed to griddata().
        '''
        if y is None:
            y = self.dye_v - self.diode_v
        tc = np.where((self.t < tl) ^ (self.t < th))
        ts = self.t[tc]
        if ngridx == None:
            ngridx = len(ts) * 10
        if ngridy == None:
            ngridy = len(y) * 10
        xi = np.linspace(np.min(ts), np.max(ts), ngridx)
        yi = np.linspace(np.min(y), np.max(t), ngridy)
        xx, yy = np.meshgrid(xi, yi)
        zg = interp.griddata((self.zc[0], self.zc[1]), self.zc[2], (xx, yy),
                **kwargs)
        plt.contuor(xi, yi, zg, 15, linewidths=0.5, colors='k')
        plt.contourf(xi, yi, zg, 15)
        return


    def specMesh(self, fl, fh, y = None,
            ngridx = 100, ngridy = 100, log = False, **kwargs):
        '''
            Displays the mesh spectrum according to the input arguments.
            Takes the logarithm of abs(g) and display it using the 'nearest'
            algorithm for mlab.griddata.

            INPUTS:
            fl      :       Frequency to start the surface at.
            fh      :       Frequency to end the surface at.
            y       :       Velocity (y axis) measure. Defaults to 
                                dye_v - diode_v
        '''
        if y is None:
            y = self.dye_v - self.diode_v
        if log == False:
            lg = np.abs(self.g)[:,np.where((self.f < fl) ^ (self.f < fh))[0]]    
        else:
            lg = np.log10(np.abs(self.g))[:,np.where((self.f < fl) ^ (self.f < fh))[0]]
        fs = self.f[np.where((self.f < fl) ^ (self.f < fh))]
        xi = np.linspace(np.min(fs), np.max(fs), ngridx)
        yi = np.linspace(np.min(y), np.max(y), ngridy)
        xx, yy = np.meshgrid(xi, yi)
        zgc = createGrid(fs, y, lg)
        zg = interp.griddata((zgc[0], zgc[1]), zgc[2], (xx, yy), **kwargs)
        plt.contour(xi, yi, zg, 15, linewidths=0.5, colors='k')
        plt.contourf(xi, yi, zg, 15)
        return(xi, yi, zg)

    def interp1dXY(self, x, y, ngridx, ngridy):
        '''
            Creates a grid of points (ngrix, ngridy) in dimensions with one
            dimensional interpolations. This is different than the typical mesh
            grid application. First, the x coordinate is interpolated for each
            row. Then the y axis is interpolated along each interpolated grid
            column.

            The y axis is generally the nonregularaly spaced one.

            This is done in this fashion since we don't need a fancy 2D qhull
            type interpolation - we only need to see the first order surface
            traced out by these data.
        '''

    def phaseMesh(self, fl, fh, y = None, ngridx = None, ngridy = None,
            **kwargs):
        '''
        Displays the mesh of phase plots.
        '''

def createGrid(x, y, z):
    '''
    Helper function to create a coordinate set for use in a griddata structure.
    INPUTS:
    X       :   x dimension.
    Y       :   y dimension.
    Z       :   z dimension. This is a 2d array of data points with
                dimensions (len(x), len(y))
    '''
    xx, yy = np.meshgrid(x, y)
    Ns = xx.shape[0] * xx.shape[1]
    zc = np.vstack([np.reshape(xx, Ns), np.reshape(yy, Ns), np.reshape(z, Ns)])
    return zc
