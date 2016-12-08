# Class file for handling data sets from f0 and f1 measurements using the lock in.
# This is to make life a little bit easier with loading / unloading things.
import lockin
import numpy as np
import clean

class xcorr:
    def __init__(self):
        fls = get_filelists()
        self.f0_0kHz_fl = fls[0]
        self.f0_20kHz_fl = fls[1]
        self.f1_20kHz_h1_fl = fls[2]
        self.f1_20kHz_h2_fl = fls[3]
        self.data_path = '/home/sean/Data/10-24-2016/'

    def createGrids(self):
        [self.f0_0kHz] = self.clean_lockin_list(self.f0_0kHz_fl)
        [self.f0_20kHz] = self.clean_lockin_list(self.f0_20kHz_fl)
        [self.f1_20kHz_h1] = self.clean_lockin_list(self.f1_20kHz_h1_fl)
        [self.f1_20kHz_h2] = self.clean_lockin_list(self.f1_20kHz_h2_fl)


    def clean_lockin_list(self, flist):
        temp = [lockin.read_lockin(self.data_path + x) for x in flist]
        wl = [x[0] for x in temp]
        re = [x[1] for x in temp]
        im = [x[2] for x in temp]
        pow = [x[3] for x in temp]
        #We need to clean and align all these data. Do this by rounding according to the wavelength value
        [wl, re] = np.vstack([clean.clean(np.around(x[0], decimals =5), x[1]) for x in zip(wl,re)])
        [wl, im] = np.vstack([clean.clean(np.around(x[0], decimals =5), x[1]) for x in zip(wl, im)])
        [wl, pow] = np.vstack([clean.clean(np.around(x[0], decimals=5), x[1]) for x in zip(wl, p)])
        return [wl, re, im, pow]



def get_filelists():
    f0_0kHz = ['f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-0p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-10p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-11p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-12p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-13p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-14p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-15p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-16p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-17p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-1p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-2p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-3p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-4p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-5p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-6p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-7p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-8p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B100Amp_sep-9p0cm.txt',
               'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_0kHz_harm-1_9p99V_B150Amp_sep-0p0cm.txt', ]
    f0_20kHz = ['f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-0p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-10p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-11p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-12p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-13p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-14p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-15p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-16p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-17p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-1p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-2p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-3p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-4p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-5p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-6p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-7p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-8p0cm.txt',
                'f0_50W_0p25mTorr_0p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-9p0cm.txt']
    f1_20kHz_h1 = ['f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-10p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-11p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-12p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-13p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-14p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-15p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-16p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-17p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-1p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-2p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-3p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-4p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-5p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-6p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-7p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-8p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-1_9p99V_B100Amp_sep-9p0cm.txt']
    f1_20kHz_h2 = ['f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-14p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-15p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-4p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-5p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-6p0cm.txt',
               'f1_50W_0p25mTorr_0p6s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-7p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-0p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-10p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-11p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-12p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-13p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-16p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-17p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-1p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-2p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-3p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-8p0cm.txt',
               'f1_50W_0p25mTorr_1p2s_3Ord_2Hz_antenna_2-4_3-5_20kHz_harm-2_9p99V_B100Amp_sep-9p0cm.txt']
    return (f0_0kHz, f0_20kHz, f1_20kHz_h1, f1_20kHz_h2)