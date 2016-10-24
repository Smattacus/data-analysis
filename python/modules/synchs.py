#
#   Some functions and objects for reading synchronization data.
#
#
#

import numpy as np
import h5py


def parseSynchFileChannel(filename, chan = 0, offset=0x8000):
    do = np.array(h5py.File(filename, 'r')['SIS3316_Channel_' + str(chan)])
    return parseSynchData(do, offset)

def parseSynchData(synch_data, offset=0x8000):
    '''
    This routine takes an array of data from the SIS3316 and returns the averaged values and timestamps
    from the raw dataset. It assumes that the default short unsigned int dataset has no raw samples and
    2 averaged samples - this works out to 10 unsigned short words per event:

    (averages, timestamps) = parseSynchdata(synch_data, offset=0x8000)

    Args:
        synch_data:     Array from the sis3316 digitizer from the h5 file. Assumes that no
         raw samples are taken and only two averaged samples are taken.
        offset   :     default=0x8000. Offset value to convert the raw short unsigned int
            to floating voltage values.

    Returns:
        averages:       nx2 array of 2 averaged samples taken by the digitizer. Rescaled by offset to give
                            double voltages.
        timestamps:     Array of timestamp values corresponding to the averaged samples.
    '''
    #Cast as an arrow just in case it hasn't already been done.
    synch_data = np.array(synch_data)
    t3 = synch_data[1::10]
    t1 = synch_data[2::10]
    t2 = synch_data[3::10]
    #bitshift the second and third chunks leftwise to create the final timestamps.
    timestamps = np.uint64(t1) + (np.uint64(t2) << 16) + (np.uint64(t3) << 32)
    #Now, take care of the data itself:
    avs1 = synch_data[8::10]
    avs2 = synch_data[9::10]
    #Subtract offset, divide by max uint value, and rescale over 5V range:
    avs1 = (avs1.astype(float) - offset) / 0xffff * 5.0
    avs2 = (avs2.astype(float) - offset) / 0xffff * 5.0
    avs = np.vstack((avs1, avs2)).T
    return avs, timestamps


#TODO: Write an interpolator function that checks the timestamps to see if a trigger
#was missed.
def interpolateTimeStamps():
    return

