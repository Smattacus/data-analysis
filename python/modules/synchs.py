#
#   Some functions and objects for reading synchronization data.
#
#
#

import numpy as np


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
    t3 = synch_data[1::10]
    t1 = synch_data[2::10]
    t2 = synch_data[3::10]
    #bitshift the second and third chunks leftwise to create the final timestamps.
    timestamps = np.uint64(t1) + (np.uint64(t2) << 16) + (np.uint64(t3) << 32)
    #Now, take care of the data itself:
    avs1 = synch_data[8::10]
    avs2 = synch_data[9::10]
    avs1 = (avs1.astype(float) - offset) / np.iinfo(np.uint16).max
    avs2 = (avs2.astype(float) - offset) / np.iinfo(np.uint16).max
    avs = np.hstack((avs1, avs2))
    return avs, timestamps

