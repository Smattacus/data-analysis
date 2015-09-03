#Module file that implements some basic Phillips - Twomey methods.

import numpy as np

def genConstB(M):
    '''
        Generates and returns an (M-1) X M first difference matrix according to 
        19.5.3 of Numerical Recipes 3rd edition.

        B = [[-1, 1, 0, ... 0],
            [0, -1, 1, 0, ..., 0],
            ...
               [0, ..., 0, 0, -1, 1]]
    '''
    if M < 2:
        print('negative matrix size detected!')
        raise RuntimeError
    fr = np.zeros((M,1))
    fr[0] = -1
    fr[1] = 1
    fr = fr.T
    B = np.vstack(np.array([np.roll(fr, i) for i in range(M-1)]))
    return B

def SimpQuad(x):
    '''
    Generates quadrature coefficients for a simpson integration scheme. Assumes
    a uniformly distributed x matrix.
    '''
    h = x[1] - x[0]
    w = np.ones((np.size(x),1))
    w[1::2] = 4
    w[2::2] = 2
    w = w * h / 3.0
    return w.T

def chisq(xdata, ydata, xk, uk, kfunc, std, quad):
    '''
    Generate and return chi square with a given kernel function kfunc and a
    fitted solution function u.
    INPUTS:
    xdata   :   x values of measured data points.
    ydata   :   Measured data points.
    xk      :   x range used for generation of kernel function
    uk      :   Fitted solution vector.
    kfunc   :   Kernel function pointer.
    std     :   Error vector, or single value (if uniform errors)
    quad    :   Quadrature weighting coefficients for use in numerical
                integration.
    '''
    fitint = np.vstack([kfunc(xk - i) * quad for i in xdata])
    fitvals = np.dot(fitint, uk)
    chisq = np.sum(ydata - fitvals)**2 / std**2
    return chisq

def LinReg(x, kernel, xdata, ydata, weight, minrule='Second', *args):
    '''
    Implements a standard Phillips - Twomey linear regularization.
    (A.T * A + w * H) f = A.T * b
    INPUTS:
    x       :   x array of the known integral kernel.
    kernel  :   Kernel function. Of form kernel(x, *args)
    xdata   :   X positions of measured data points.
    ydata    :   Measured data points.
    weight  :   Weight value.
    *args   :   Additional arguments passed to integral kernel.
    '''
    #Use the xdata points as the offsets for generating the design matrix.
    A = genDesign(x, xdata, kernel, *args)
    if minrule=='Second':
        B = genLinB(np.size(x))
    H = np.dot(B.T, B)
    blm = np.dot(A.T, A) + w * H
    brm = np.dot(A.T, ydata)
    sol = la.lu_solve(la.lu_factor(blm), brm)
    return sol
    

def genDesign(x, offsets, func, quad='Simpson', *args):
    '''
    Generate a design matrix A for a Phillips - Twomey linear regularization.
    Used to set up the linear equation for the inverse problem:
    (A.T * A + l * H ) f = A.T * b
    INPUTS:
    x       :   x array of data points.
    offsets :   Array of x values that are the measured data point positions.
    func    :   Function of the known underlying physical process.
    quad    :   Any given quadrature coefficient array for the dataset. Simpson
    is usually a good choice.
    '''
    if quad=='Simpson':
        q = SimpQuad(x)
    else:
        print('Implement other quadrature methods!')
        raise
    A = np.vstack([func(x - i, *args) * q for i in offsets])
    return A

def genLinB(M):
    '''
    Generates and returns an (M-2) x M second difference matrix according to
    19.5.11 of Numerical Recipes 3rd edition.

    So B = [[-1, 2, -1, 0, ..., 0],
            [0, -1, 2, -1, 0, ..., 0],
            ...
            [0, ..., 0, -1, 2, -1]]
    '''
    if M < 3:
        print('Negative matrix size detected!')
        raise RuntimeError
    fr = np.zeros((M,1))
    fr[0] = -1
    fr[1] = 2
    fr[2] = -1
    fr = fr.T
    B = np.vstack(np.array([np.roll(fr, i) for i in range(M-2)]))
    return B

def genQuadB(M):
    '''
    Generates and returns an (M-3) x M difference matrix according to 19.5.14
    of Numerical Recipes 3rd edition.

    B = [[-1, 3, -3, 1, 0, ..., 0],
        [0, -1, 3, -3, 1, 0, ..., 0],
        ...
        [0, ..., 0, -1, 3, -3, 1]]
    '''
    if M < 4:
        print('Negative matrix size detected!')
        raise RuntimeError
    fr = np.zeros((M,1))
    fr[0] = -1
    fr[1] = 3
    fr[2] = -3
    fr[3] = 1
    fr = fr.T
    B = np.vstack(np.array([np.roll(fr, i) for i in range(M-3)]))
    return B

