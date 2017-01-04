#!/usr/bin/python2

#This is a quick script to generate the .m and .sh files to run MATLAB jobs,
#given a listing of files.

import sroh

dprefix = '/Users/smattingly/Data/'
xcprefix = '/Users/smattingly/XCMeans/Cuts/'

#Template for the sh file.

shtemplate = """#!/bin/sh
#Script to generate XCMeans and spec byproducts with cutting.
#Data from DATE

#$ -pe smp 12
#$ -N XCMC_DATE_DIODE_DYE
#$ -o XCMC_DATE_DYE.log -q UI -cwd
#$ -M sean-mattingly@uiowa.edu
#$ -m abe

matlab -nodesktop -nosplash -r 'saveXCM_MDATE_DIODE_DYE'"""

#Template for the matlab .m file.

mtemp = """%Script to print phases to a file in MATLAB.
%Template
addpath('~/Programs/data-analysis/MATLAB/');
xcmpath = 'XCM_PREFIX/DATE/DIODE/';
dyename = 'DYE';
datapath = sprintf('DATAPATH/DATE/DIODE/DYE/', dyename);
savename = sprintf('%s/%s.mat', xcmpath, dyename);
"""

cuttemp = """saveXCMeanBoxFindPhase(datapath, savename, 1e5, 1e6, CUTVAL);
exit;"""

nctemp = """saveXCMeanBoxFindPhase(datapath, savename, 1e5, 1e6);
exit;"""

fcl = open('/home/sean/data/SeanRunOne/XCMeanScripts/file_cutlist.txt', 'r')
shs = []
shtitles = []
ms = []
mtitles = []
dcs = []
sht = 'save_XCM_DATE_DIODE_DYE.sh'
mt = 'saveXCM_MDATE_DIODE_DYE.m'
for x in fcl:
    dcs.append(x.split())

#We now have a listing of data files and directories with cut values.
for x in dcs:
    date = x[0].split('/')[2]
    mdate = date.replace('-', '_')
    diode = x[0].split('/')[3]
    dye = x[0].split('/')[4].split('.')[0]
    shcurr = shtemplate.replace('MDATE', mdate).replace('DATE', date).replace('DYE',
    dye).replace('DIODE', diode)
    mcurr = mtemp.replace('XCM_PREFIX', xcprefix).replace('DATE',
            date).replace('DYE', dye).replace('DATAPATH',
            dprefix).replace('DIODE', diode)
    if len(x) > 1:
        mcurr = mcurr + cuttemp.replace('CUTVAL', x[1])
    else:
        mcurr = mcurr + nctemp
    shs.append(shcurr)
    ms.append(mcurr)
    shtitles.append(sht.replace('DIODE', diode).replace('DYE',
        dye).replace('DATE', date))
    mtitles.append(mt.replace('DIODE', diode).replace('DYE',
        dye).replace('MDATE', mdate).replace('DATE', date))

#Now we can just zip these four lists togother and write all these files.
for x in zip(shs, ms, shtitles, mtitles):
    date = x[2].split('_')[2]
    fsh = open(date + '/' + x[2], 'w')
    fm = open(date + '/' + x[3], 'w')
    fsh.write(x[0])
    fm.write(x[1])
    fsh.close()
    fm.close()
