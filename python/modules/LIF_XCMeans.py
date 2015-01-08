#This is a file containing major functions from MATLAB used for the purpose of
#reducing data from the Skiff lab HLMX experiment.

#This is to be used to replace the MATLAB stuff finally.

def genAllTPBoxFindPhase(path, Fc, Fa):
    '''
Function to generate ALL the downsampled TOP and BOT arrays of LIF + Noise
and noise data. 

[TOP1, BOT1, TOP2, BOT2] = genAllTPBoxFindPhase(path, Fc, Fa)

This function assumes that the same phase is to be used across all files
in the directory pointed to by path. It generates a straight square wave
with phase and frequency given by phase and freq, respectively. This
square wave is used to demodulate the data.

This function finds all available files by searching for .h5 files. It
will use all of them in the directory, so delete or move unwanted ones.

INPUTS:
   path - Path to a directory of h5 files.
   Fc - Frequency of the square wave function. Generally, this is the
            laser chop frequency
   Fa   - Acquisition frequency.

OUTPUTS:
   TOP1 - NxM array of LIF +plasma data. N = number of files, M = # of points, PMT 1.
   BOT1 - NxM array of plasma data, PMT 1.
   TOP2 - NxM array of LIF + plasma data, PMT 2.
   BOT2 - NxM array of plasma data, PMT 2.
    '''


