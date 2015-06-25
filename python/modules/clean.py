#function for cleaning an x and y set of data that isn't regularly spaced
#on the x axis.

#Based on a program written by F W Skiff in MATLAB.
import numpy as np

def clean(x, y):
    """
        Regularizes the distances in x of the data. Interpolates gaps in the y
        data or averages repeating values.

        This will sort the input data - it doesn't need to be presorted.
    """
    xi= np.argsort(x)
    x = x[xi]
    y = y[xi]
    (u, ui) = np.unique(x, return_index=True)
    di = np.diff(ui)
    #Average redundant data.
    L0 = [np.mean(y[ui[i]:ui[i]+di[i]]) for i in range(np.size(di))]  
    L0 = L0 + [np.mean(y[ui[np.size(di):]])]
    L0 = np.array(L0)
    #Interpolate
    xls = np.linspace(u[0], u[-1], np.size(u))
    L0i = np.interp(xls, u, L0)
    return (xls, L0i) 
