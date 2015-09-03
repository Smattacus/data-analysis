#Module of various function for streamlining loading of data.
#Contains function get arrays and addresses for data files.
import glob
import numpy as np

def getDiodePoints():
    #Hand checking, these match the entries in my logbook.
    diode_points = [
    668.61091,  
    668.61095,  
    668.61100,  
    668.61105,  
    668.61110,  
    668.61117,  
    668.61114,  
    668.61119,  
    668.61127,  
    668.61131,  
    668.6114024, #668.6114024 is from notebook. Rounding for numerical reasons.
    668.61142,  
    668.61147,  
    668.61146,  
    668.61150, 
    668.61167, 
    668.61179, 
    668.61171,
    668.61180,
    668.61199, 
    668.61185,
    ]
    return diode_points

#Rels
#From the dye handconvolve fits
def getRels():
    #These shifts are given relative to the first day's fitted Gaussian.
    #i.e., (fits) - (first gaussian)
    rels = np.array([  0.00000000e+00, #4-28-15 [0]
             5.32548738e-05, #4-23-15 [1]
            -2.78384236e-04, #4-17-15 - higher P (this was used for data) [2]
            -3.29720871e-04, #4-17-15 - lower P (this was not used) [3]
            -1.50167993e-04, #4-14-15 [4]
            -4.55705386e-04, #4-10-15 [5]
            -4.12182305e-04, #4-08-15 [6]
            -5.81653718e-04, #4-06-15 [7]
            -3.04319509e-04, #4-03-15 [8]
            -4.49631641e-04, #4-07-15 [9]
            -4.42360542e-04, #4-09-15 [10]
            -3.59077275e-04, #4-13-15 [11]
            -1.72130652e-04, #4-15-15 [12]
            -8.51021464e-05, #4-21-15 [13]
             7.57091323e-06, #4-22-15 [14]
             4.48569494e-05, #4-24-15 [15]
            -7.69731422e-05]) #4-27-15 [16]
    return rels

def getDyeIVDFFiles():
    ivdflist= [ 
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr28_2015.txt', # -10
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr23_2015.txt.txt', # -9, -8
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr17_2015_higherP.txt', # -7, -6
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr14_2015.txt', # +5, -5
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr10_2015.txt', # -3, -4
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr8_2015.txt', # -2
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr6_2015.txt', # -1
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr3_2015.txt', # 0
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr7_2015.txt', # +1, +2
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr9_2015.txt', # +3
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr13_2015.txt', # +4
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr15_2015.txt', # +6
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr21_2015.txt', # +7
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr22_2015.txt', # +8
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr24_2015.txt', # +9
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DYE_Apr27_2015.txt' # +10
            ]
    return ivdflist

def getDiIVDFFiles():
    #Unfortunately, the data file for April 8th is gone (overwritten accidentally).
    dilist = [ '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr28_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr23_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr17_2015_higherP.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr17_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr14_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr10_2015.txt',
                 '',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr6_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr3_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr7_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr9_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr13_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr15_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr21_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr22_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr24_2015.txt',
            '/home/sean/data/SeanRunOne/XCMeans/IVDFS/DIODE_Apr27_2015.txt' ]
    return dilist

def getDataDirs():
    dirlist = [ '/home/sean/data/SeanRunOne/XCMeans/04-28-15/Diode_N10Step',
    '/home/sean/data/SeanRunOne/XCMeans/04-23-15/Diode_N9Step',
    '/home/sean/data/SeanRunOne/XCMeans/04-23-15/Diode_N8Step',
    '/home/sean/data/SeanRunOne/XCMeans/04-17-15/Diode_N7Step',
    '/home/sean/data/SeanRunOne/XCMeans/04-17-15/Diode_N6Step',
    '/home/sean/data/SeanRunOne/XCMeans/04-14-15/Diode_N5Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-10-15/Diode_N4Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-10-15/Diode_N3Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-08-15/Diode_N2Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-06-15/Diode_N1Step',
    '/home/sean/data/SeanRunOne/XCMeans/04-03-15', 
    '/home/sean/data/SeanRunOne/XCMeans/04-07-15/Diode_P1Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-07-15/Diode_P2Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-09-15/Diode_P3Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-13-15/Diode_P4Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-14-15/Diode_P5Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-15-15/Diode_P6Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-21-15/Diode_P7Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-22-15/Diode_P8Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-24-15/Diode_P9Step', 
    '/home/sean/data/SeanRunOne/XCMeans/04-27-15/Diode_P10Step']
    return dirlist

def getShiftIndices():
    shift_indices = (np.array([0, 1, 1,
    2, 2, 4, 5, 5, 6, 7, 8, 9, 9, 10, 11, 4, 12, 13, 14, 15, 16]),)
    return shift_indices 

def getDiodeWavelengths():
    di_p = getDiodePoints()
    s_i = getShiftIndices()
    shift = getRels()[s_i]
    new_di_p = [x[0] + x[1] for x in zip(di_p, shift)]
    return new_di_p

def getDyeWavelengths(): 
    dl = getDataDirs()
    shift_indices = getShiftIndices()
    shifts = getRels()[shift_indices]
    dye_wls = []
    for x in dl:
        temp = glob.glob(x + '/Dye_*.mat')
        temp.sort()
        #split ops: split('/') for filename, '_' for 611p?????.mat, '.' for 611p?????,
        #replace('p', '.').
        dye_wls.append([float(y.split('/')[-1].split('_')[-1].split('.')[0].replace('p','.'))
            for y in temp])
    return dye_wls
