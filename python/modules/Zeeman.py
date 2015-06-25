#Functions for calculating the splitting of Zeeman lines.
#Assumes anomolous Zeeman effect.

h_cgs = 6.62606957e-27
c_cgs = 2.99792456e10
o_cgs = 9.27400968e-21

import numpy as np
from numpy import linspace
import clean
import scipy as sp
import spec

#NOTE this value is for adjusting the windowing in the functions
#ZeemanFitFunc_VaryB and ZeemanFitFunc_VaryV. Adjust in python environment if
#needed.
#BE CAREFUL.
freq_end = 1000
fnoise_start = 1500
fnoise_end = 7000
print('Using global variable: freq_end . Be careful.')

RS_Ls = {'S':0.0, 'P':1.0, 'D':2.0, 'F':3.0, 'G':4.0, 'H':5.0, 'I':6.0,
'K':7.0,'L':8.0,'M':9.0}

def RS(S, L, J):
    """
    Defines a Russel - Saunders dict for the three Russel Saunders terms.
    
    Just enter an RS term exactly as it's written.
    INPUTS:
        S   :   Spin multiplicity.
        L   :   Orbital quantum number in spectroscopy notation. Enter
        capitalized orbital (S, P, D, F, ...)
        J   :   Total angular momentum.
    """
    return {'S': np.double(S), 'L':np.double(RS_Ls[L]), 'J':np.double(J)}

def Splittings(E_o, RSterm, B):
    """
        Returns a list of the energy levels of the Zeeman split line. This
        assumes a weak field Zeeman splitting with Russel - Saunders coupling.

        Splittings(E_o, RSterm, B)

        INPUTS:
            E_o     :   Initial Energy level
            RSTerm  :   Dict containing the Russel - Saunders coupling terms.
            B       :   Magnetic field in Gauss.
        OUTPUTS:
            Esplits :   (numeric array) Energy levels of split lines
            Msplits :   (numeric array) Angular momentum values of split lines.

        This calculation is done in CGS. 
    """
    #Use the Russel - Saunders term to get S, L, and J.
    S = (RSterm['S'] - 1) / 2
    L = RSterm['L']
    J = RSterm['J']
    Mz = linspace(-J, J, num=(2 * J + 1))
    Ez = []
    g = 1 + (J * (J+1) + S*(S+1) - L*(L+1)) / (2 * J * (J+1))
    #print(g)
    #print(o_cgs * B * 0.5 * g / (h_cgs * c_cgs))
    for i in Mz:
        Ez.append(o_cgs * B * i * g / (h_cgs * c_cgs) + E_o)
    return (np.array(Ez), Mz)

def Transitions(Ez_i, Mz_i, RS_i, Ez_f, Mz_f, RS_f):
    """
    Returns arrays of the allowed transitions in energy (cm-1), wavelength
    (nm), and polarization (+1, -1, or 0).

    Takes as inputs two lines from the Zeeman.Splittings routine. It will check
    selection rules before running.

    INPUTS:
        Ez_i    :   (array, numeric) Energy levels of the initial state.
        Mz_i    :   (array, numeric) space quantization values of the initial
                        energy level.
        RS_i    :   (dict) Russel - Saunders term for the initial level.
        Ez_f    :   (array, numeric) Energy levels of the final state.
        Mz_f    :   (array, numeric) space quantization values of the final
                        energy level.
        RS_i    :   (dict) Russel - Saunders term for the final level.
    OUTPUTS:
        Etrans  :   (numeric) Array of transition energies in cm-1.
        Wtrans  :   (numeric) Array of wavelengths (nm)
        Polarization    :   (numeric) Array of polarization values.
                            +1: right hand circular polarize (sigma+)
                            -1: Left hand circular polarize (sigma-)
                            0:  Linear (pi).
        Mpairs  :   (list) List of tuples corresponding to the initial and
                    final spatial quantization levels. (Mi, Mf)

    """
    Si = (RS_i['S'] -1) / 2
    Li = RS_i['L']
    Ji = RS_i['J']
    Sf = (RS_f['S'] - 1) / 2
    Lf = RS_f['L']
    Jf = RS_f['J']
    Js = Jf - Ji
    if Js != 1 and Js != 0 and Js != -1:
        print('Selection rule for J in RS term violated! Double check your inputs!')
        return
    #Handle the transitions over three separate loops: dm = 0, 1, -1.
    if Js == 0:
        #delta_m = 0 is not allowed by selection rules if Js = 0.
        Etrans_0 = np.array([])
        Wtrans_0 = np.array([])
        Ptrans_0 = np.array([])
    else:
        #Do delta_m = 0
        if np.size(Mz_i) > np.size(Mz_f):
            #delta_m = 0 is constrained to Mz_f m values.
            sdiff = np.size(Mz_i) - np.size(Mz_f)
            Etrans_0 = Ez_i[sdiff / 2:-sdiff/2] - Ez_f
            Mpairs_0 = zip(Mz_i[sdiff / 2: -sdiff/ 2], Mz_f)
        elif np.size(Mz_i) < np.size(Mz_f):
            #delta_m = 0 is constrained to Mz_o m values.
            sdiff = np.size(Mz_f) - np.size(Mz_i)
            Etrans_0 = Ez_i - Ez_f[sdiff/2:-sdiff/2]
            Mpairs_0 = zip(Mz_i, Mz_f[sdiff/2:-sdiff/2]) 
        else:
            #They are the same size.
            Etrans_0 = Ez_i - Ez_f
            Mpairs_0 = zip(Mz_i, Mz_f)
        Wtrans_0 = 1 / Etrans_0
        Ptrans_0 = Wtrans_0 * 0
    #Handle the next transitions.
    #delta_m = +1.
    Etrans_p1 = np.array([])
    Wtrans_p1 = np.array([])
    Ptrans_p1 = np.array([])
    Mpairs_p1 = []
    for i in (Mz_i + 1):
        if i in Mz_f:
            #The transition has an allowed final result: calculate.
            mi = i - 1
            mf = i
            Et = Ez_i[np.where(Mz_i == mi)] - Ez_f[np.where(Mz_f == mf)] 
            Etrans_p1 = np.append(Etrans_p1, Et)
            Wtrans_p1 = np.append(Wtrans_p1, 1 / Et)
            Ptrans_p1 = np.append(Ptrans_p1, 1)
            Mpairs_p1 = Mpairs_p1 + [(mi, mf)]
    #delta_m = -1
    Etrans_n1 = np.array([])
    Wtrans_n1 = np.array([])
    Ptrans_n1 = np.array([])
    Mpairs_n1 = []
    for i in (Mz_i - 1):
        if i in Mz_f:
            mi = i + 1
            mf = i
            Et = Ez_i[np.where(Mz_i == mi)] - Ez_f[np.where(Mz_f == mf)]
            Etrans_n1 = np.append(Etrans_n1, Et)
            Wtrans_n1 = np.append(Wtrans_n1, 1 / Et)
            Ptrans_n1 = np.append(Ptrans_n1, -1)
            Mpairs_n1 = Mpairs_n1 + [(mi, mf)]
    Etot = np.append([Etrans_0, Etrans_p1], Etrans_n1)
    Wtot = np.append([Wtrans_0, Wtrans_p1], Wtrans_n1)
    Ptot = np.append([Ptrans_0, Ptrans_p1], Ptrans_n1)
    Mtot = Mpairs_0 + Mpairs_p1 + Mpairs_n1
    return (Etot, Wtot, Ptot, Mtot)

def Amplitudes(RSi, RSf, Mt):
    """
        Calculates the amplitudes of the transitions given by
        Zeeman.Transitions()
        
        Inputs:
            RSi     :   Russel Saunders term for initial level.
            RSf     :   Russel Saunders term for final level.
            Mt      :   List of tuples containg initial and final
                        space quantization values. Use output from
                        Z.Transition().
        Outputs:
            Amps    :   (array) Numeric Array of relative intensities of the
            transitions.
    """
    Jt = RSf['J'] - RSi['J']
    Ji = RSi['J']
    Jf = RSf['J']
    Amps = []
    if Jt == 0:
        for i in Mt:
            if (i[1] - i[0]) == 1:
                Amps = Amps + [(Ji - i[0]) * (Ji + i[0] + 1)]
            elif (i[1] - i[0]) == 0:
                Amps = Amps + [2 * (i[0] ** 2)]
            elif (i[1] - i[0]) == -1:
                Amps = Amps + [(Ji + i[0]) * (Ji - i[0] + 1)]
            else:
                print('Selection rule violated! J == 0')
                raise
                return -1
    elif Jt == 1:
        for i in Mt:
            md = i[1] - i[0]
            mi = i[0]
            if md == 1:
                Amps = Amps + [(Ji + mi + 1) * (Ji + mi + 2)]
            elif md == 0:
                Amps = Amps + [2 * ((Ji + 1)**2 - mi**2)]
            elif md == -1:
                Amps = Amps + [(Ji - mi + 1) * (Ji - mi + 2)]
            else:
                print('Selection rule violated! J == 1')
                raise
    elif Jt == -1:
        for i in Mt:
            md = i[1] - i[0]
            mi = i[0]
            if md == 1:
                Amps = Amps + [(Ji - mi) * (Ji - mi - 1)]
            elif md == 0:
                Amps = Amps + [ 2 * (Ji**2 - mi**2)]
            elif md == -1:
                Amps = Amps + [(Ji + mi) * (Ji + mi - 1)]
            else:
                print('Selection rule violated! J == -1')
                raise
    else:
        print('Selection rule violated for J! returning.')
        raise
    return Amps

def Linewidths_Diode(wavelengths, P, velocity=0):
    """
        Returns an array corresponding to the input wavelength values of the
        sum of all the Lorentzian functions for the Ar II transition
        corresponding to the Diode transition:

        4F_(7/2) <- 4D_(5/2)^0

        With an input laser wavelength of 668.61379.

        This is done for the Skiff plasma experiment, assuming a Zeeman
        splitting with a 1kG electric field.

        INPUTS:

            wavelengths :   (numeric) Array of wavelengths to evaluate the
                            Lorentzian functions at. Enter in nanometers.
            P           :   (int)   Integer denoting the polarization to
                            generate the spectrum for.
            velocity=0  :   Streaming flow. cm/s with positive towards the
                            observer. Optional.

        OUTPUTS:
            spectrum    :   (numeric) Sum of all the lorentzians for the given
                            input wavelength values.

        The S(lambda) output is of the form:

        S(lambda) = A_1 * L_1(lambda) + A2 * L2(lambda) + ...

        Where A_i corresponds to the relative amplitude of the Zeeman
        transition while L_i is the natural linewidth.
    """
    RSi = RS(4, 'D', 5 * 0.5)
    RSf = RS(4, 'F', 7 * 0.5)
    Ei = 157673.4134
    Ef = 142717.0965

    ti = 1.0/1.326e8

    (Ez_f, Mz_f) = Splittings(Ef, RSf, 1e3)
    (Ez_i, Mz_i) = Splittings(Ei, RSi, 1e3)

    (Et, Wt, Pt, Mt) = Transitions(Ez_i, Mz_i, RSi, Ez_f, Mz_f, RSf)
    At = Amplitudes(RSi, RSf, Mt)
    
    Zinfo = np.array([Et, Wt * 1e7, At, Pt, np.array(Mt)[:,0], np.array(Mt)[:,1]]).transpose()

    #Zinfo_P now contains the Energy | Wavelength (nm) | Amplitude |
    #Polarization | Mz_I | Mz_F 
    #of the transitions with Polarization = P.
    Zinfo_P = Zinfo[np.where(Zinfo[:,3] == P), :][0]
    Amps = Zinfo_P[:,2]
    Cents = Zinfo_P[:,1]

    #Amplitudes and centers are all we need.
    Ltot = np.zeros(np.size(wavelengths))
    for x in zip(Amps, Cents):
        Ltot += ZL(x[0], x[1], (ti,), wavelengths, velocity = velocity)
    return Ltot

def ZL(Amp, Cent, ts, wavelengths, velocity=0):
    """
        Helper function for generating a Lorentzian for the natural linewidth.
        INPUTS:
            Amp         :   Relative amplitude of the line.
            Cent        :   Center of the line in nm.
            ts          :   Tuple containing the mean lifetimes of the final and
                        initial energy levels. Only one element in the tuple is fine.
            wavelengths :   Wavelengths to return the evaluted lorentzian over.
                            (nm)
            velocity=0  :   (float) (Optional) Add a streaming flow in cm/s to
                            result. This shifts all lines by a flat amount
                            given by
                            sigma = sigma_0 * ( 1 + u / c)
                            where u is the velocity of the particle moving
                            towards the observer, and sigma and sigma_0 are in
                            cm-1 units.
        OUTPUTS:
            L           :   Lorentzian function.

        NOTES:
        The generating Lorentzian is according to Condon & Shortley, with 

        L(lambda) = 1 / pi * (1 / sigma_0 / (1 + (sigma/sigma_0^2)))

        and

        sigma_0 = 1 / (4 pi c) (1 / tau(A) + 1 / tau(B))

        Where tau(A) and tau(C) are the mean lifetimes of the initial and final
        states. If there are not transitions from the lower state (i.e. it's a
        ground or stable state) then it can be safely dropped.
    """
    wlcm = 1/(wavelengths * 1e-7)
    wlcm = wlcm * (1 + velocity / c_cgs)
    IC = 1 / (Cent * 1e-7)
    s0 = 0
    for x in ts:
        s0 += 1 / (4 * np.pi * c_cgs) * 1 / x
    S = Amp / np.pi * (1 / s0) / (1 + ((IC - wlcm) / s0)**2)
    return S

def Zeeman_Lines(RSi, RSf, Ei, Ef, wavelengths, P, B, tls, velocity = 0):
    """
        Returns the broadened Zeeman lines corresponding to a certain spectrum.

        INPUTS:
        RSi     :   RS term for the initial level.
        RSf     :   RS term for the final level.
        Ei      :   Energy (cm-1) for the initial level.
        Ef      :   Energy (cm-1) for the final level.
        wavelengths :   Wavelength (nm) for the levels.
        P       :   Polarization of desired lines. 0, 1, -1.
        tls     :   tuple of the mean lifetime of the initial and final levels
                    (ti, tf). Just one element is fine (ti or tf).
    """
    (Ez_f, Mz_f) = Splittings(Ef, RSf, B)
    (Ez_i, Mz_i) = Splittings(Ei, RSi, B)
    (Et, Wt, Pt, Mt) = Transitions(Ez_i, Mz_i, RSi, Ez_f, Mz_f, RSf)
    At = Amplitudes(RSi, RSf, Mt)
    Zinfo = np.array([Et, Wt * 1e7, At, Pt, np.array(Mt)[:,0], np.array(Mt)[:,1]]).transpose()
    Zinfo_P = Zinfo[np.where(Zinfo[:,3] == P), :][0]
    Amps = Zinfo_P[:,2]
    Cents = Zinfo_P[:,1]

    #Amplitudes and centers are all we need
    Ltot = np.zeros(np.size(wavelengths))
    for x in zip(Amps, Cents):
        Ltot += ZL(x[0], x[1], tls, wavelengths, velocity=velocity)
    return Ltot

def ZeemanFitFunc_VaryB(RSi, RSf, Ei, Ef, P, B, tls, wavelengths, s):
    """
        Function which takes a Zeeman spectrum, deconvolves it with the natural
        Zeeman linewidth according to B and P, and returns the deconvolved
        spectrum. This is to be used with an optimization problem.

        See also: ZeemanFit_VaryV, which postulates that there is a streaming
        flow in the plasma.
        
        INPUTS:
        RSi     -   Russel Saunders term for the initial energy level.
        RSf     -   Russel Saunders term for the final energy level.
        Ei      -   Energy for the initial energy level.
        Ef      -   Energy for the final energy level.
        P       -   Polarization of Zeeman lines to fit (-1, 0, 1).
        B0      -   Initial guess for magnetic field in Gauss.
        tls     -   Tuple of (ti, tf) for the mean lifetimes of the starting
                    and ending state. Just one entry in the tuple is fine.
        wavelengths -   Wvelengths of data.
        spec    -   Measured LIF spectrum.
        OUTPUTS
        x       -   Cleaned wavelengths.
        yn      -   Cleaned and normalized spectrum.
        lns_clnn-   Cleaned and normalized natural lines.
        ivdf    -   Resulting deconvolved spectrum.
    """
    #Try deconvolving using the usual process.
    #Normalize:
    (x, y) = clean.clean(wavelengths, s)
    yn = y / np.sqrt(np.dot(y, y))
    #Pad the Zeeman array so that the Zeeman lines aren't
    #truncated by numerical errors.
    #Do this by adding high frequency elements to the FFT
    #of the original data array.
    pad_length = 5000
    [f, gs] = spec.spec(yn, x[1] - x[0])
    gs = np.pad(gs, pad_length, mode='constant', constant_values=0)
    xp = linspace(x[0], x[-1], np.size(gs))
    lns_cln = Zeeman_Lines(RSi, RSf, Ei, Ef, xp, P, B, tls)
    lns_clnn = lns_cln / np.sqrt(np.dot(lns_cln, lns_cln))
    lns_clnn = lns_clnn * np.max(yn) / np.max(lns_clnn)
    [junk, gl] = spec.spec(lns_clnn, x[1] - x[0])
    #remove the HF padding.
    gs = gs[pad_length:-pad_length]
    gl = gl[pad_length:-pad_length]
    [junk, lns_clnn] = spec.ispec(gl, f[1] - f[0])
    #fi = find outside of window.
    fi = np.where(np.abs(f) > freq_end)
    #fni = noise estimation interval
    fni = np.where((np.abs(f) > fnoise_start) - (np.abs(f) > fnoise_end))
    #fwi = indices of window to keep.
    fwi = np.where(np.abs(f) <= freq_end)
    gdc = np.abs(gs) / np.abs(gl)
    gdc[fi] = 0
    [trash, ivdf] = spec.ispec(gdc, f[1] - f[0])
    ivdf = sp.fftpack.fftshift(ivdf)
    #Get the errors too.
    err = np.mean(np.abs(gs[fi]))
    sig_real = np.sqrt(err**2 / np.size(fwi) * np.sum(1 / gl[fwi]**2) *
        np.size(fwi) / np.size(f))
    errs = np.abs(sig_real)
    return (x, yn, lns_clnn, ivdf, errs)

def ZeemanFitFunc_VaryV(RSi, RSf, Ei, Ef, P, V, tls, wavelengths, s):
    """
        Function which takes a Zeeman spectrum, deconvolves it with the natural
        Zeeman linewidth according to B and P, and returns the deconvolved
        spectrum. This is to be used with an optimization problem.

        See also: ZeemanFit_VaryB, which postulates that there is a magnetic
        field other than 1kG in the plasma.
        
        INPUTS:
        RSi     -   Russel Saunders term for the initial energy level.
        RSf     -   Russel Saunders term for the final energy level.
        Ei      -   Energy for the initial energy level.
        Ef      -   Energy for the final energy level.
        P       -   Polarization of Zeeman lines to fit (-1, 0, 1).
        V       -   Streaming flow. cm/s defined to be positive towards the
                    observer (laser).
        tls     -   Tuple of (ti, tf) for the mean lifetimes of the starting
                    and ending state. Just one entry in the tuple is fine.
        wavelengths -   Wvelengths of data.
        spec    -   Measured LIF spectrum.
        OUTPUTS
        x       -   Cleaned wavelengths.
        yn      -   Cleaned and normalized spectrum.
        lns_clnn-   Cleaned and normalized natural lines.
        ivdf    -   Resulting deconvolved spectrum.
    """
    #Try deconvolving using the usual process.
    #Normalize:
    (x, y) = clean.clean(wavelengths, s)
    yn = y / np.sqrt(np.dot(y, y))
    B = 1e3
    lns_cln = Zeeman_Lines(RSi, RSf, Ei, Ef, x, P, B, tls, velocity=V)
    lns_clnn = lns_cln / np.sqrt(np.dot(lns_cln, lns_cln))
    [f, gl] = spec.spec(lns_clnn, x[1] - x[0])
    [f, gs] = spec.spec(yn, x[1] - x[0])
    fi = np.where(np.abs(f) > freq_end)
    fni = np.where((np.abs(f) > fnoise_start) - (np.abs(f) > fnoise_end))
    fwi = np.where(np.abs(f) <= freq_end)
    gdc = np.abs(gs) / np.abs(gl)
    gdc[fi] = 0
    [trash, ivdf] = spec.ispec(gdc, f[1] - f[0])
    ivdf = sp.fftpack.fftshift(ivdf)
    #Get the errors too.
    err = np.mean(np.abs(gs[fi]))
    sig_real = np.sqrt(err**2 / np.size(fwi) * np.sum(1 / gl[fwi]**2) *
        np.size(fwi) / np.size(f))
    errs = np.abs(sig_real)
    return (x, yn, lns_clnn, ivdf, errs)

def ZeemanFit_VaryB_Diode(fitvals, wavelengths, s):
    """
    This is the function for fitting the Diode to pass to the
    scipy.optimize.minimize function.

    INPUTS:
    fitvals :   List of (B, mu, sigma, amp) in order to minimize the function
                resulting from the deconvolved LIF spectrum.
                B - Magnetic field
                mu - Mean of fitted gaussian
                sigma - width of fitted gaussian
                amp - amplitude of fitted gaussian.
    wavelengths :   Wavelength values. Will be cleaned.
    s       :   Measured LIF spectrum.
    Fits the Gaussian

    A * exp (-(x - mu)**2 / (2 * sigma**2))

    and returns the sum of squares of (np.abs(ivdf) - amp * np.exp(-(x - mu)**2 / (2 * sigma**2)))

    for minimization purposes.
    """
    B = fitvals[0]
    mu = fitvals[1]
    sigma = fitvals[2]
    amp = fitvals[3]
    C = fitvals[4]
    P = -1
    RSi = RS(4, 'D', 5 * 0.5)
    RSf = RS(4, 'F', 7 * 0.5)
    Ei = 157673.4134
    Ef = 142717.0965
    ti = 1.0/1.326e8
    (x, yn, lns_clnn, ivdf, errs) = ZeemanFitFunc_VaryB(RSi, RSf, Ei, Ef, P, B, (ti,), wavelengths, s)
    ivdf = np.abs(ivdf)
    gmin = (ivdf - (amp * np.exp(-(x - mu)**2 / (2 * sigma**2)) + C)) / errs
    return np.sum(gmin**2)

def ZeemanFit_VaryV_Diode(fitvals, wavelengths, s):
    """
    Routine to pass to opt.minimize for adding a constant flow to the plasma.

    Here, the sign convention is taken to be positive for the flow moving
    towards the laser (cf Condon & Shortley p.139.

    sigma = sigma_0 * (1 + u / c)

    INCOMPLETE

    """
    V = fitvals[0]
    mu = fitvals[1]
    sigma = fitvals[2]
    amp = fitvals[3]
    C = fitvals[4]
    P = -1
    RSi = RS(4, 'D', 5 * 0.5)
    RSf = RS(4, 'F', 7 * 0.5)
    Ei = 157673.4134
    Ef = 142717.0965
    ti = 1.0/1.326e8
    (x, yn, lns_clnn, ivdf, errs) = ZeemanFitFunc_VaryV(RSi, RSf, Ei, Ef, P, V, (ti,), wavelengths, s)
    ivdf = np.abs(ivdf)
    gmin = (ivdf - (amp * np.exp(-(x - mu)**2 / (2 * sigma*2)) + C)) / errs
    return np.sum(gmin**2)

def ZeemanFit_VaryB_Dye(fitvals, wavelengths, s):
    """
    This is the function for fitting the Dye to pass to the
    scipy.optimize.minimize function.

    INPUTS:
    fitvals :   List of (B, mu, sigma, amp) in order to minimize the function
                resulting from the deconvolved LIF spectrum.
                B - Magnetic field
                mu - Mean of fitted gaussian
                sigma - width of fitted gaussian
                amp - amplitude of fitted gaussian.
    wavelengths :   Wavelength values. Will be cleaned.
    s       :   Measured LIF spectrum.
    Fits the Gaussian

    A * exp (-(x - mu)**2 / (2 * sigma**2)) + C

    and returns the sum of squares of 
    (np.abs(ivdf) - (amp * np.exp(-(x - mu)**2 / (2 * sigma**2)) + C) / err

    for minimization purposes.
    """
    B = fitvals[0]
    mu = fitvals[1]
    sigma = fitvals[2]
    amp = fitvals[3]
    C = fitvals[4]
    P = -1
    RSi = RS(2, 'F', 7 * 0.5)
    RSf = RS(2, 'G', 9 * 0.5)
    Ef = 154181.4942
    Ei = 170530.4040
    ti = 1.0/1.326e8
    (x, yn, lns_clnn, ivdf, errs) = ZeemanFitFunc_VaryB(RSi, RSf, Ei, Ef, P, B, (ti,), wavelengths, s)
    ivdf = np.abs(ivdf)
    gmin = (ivdf - (amp * np.exp(-(x - mu)**2 / (2 * sigma**2)) + C)) / errs
    return np.sum(gmin**2)

def Dye_Values():
    """
    Returns the RS terms, energy levels, and lifetimes of the standard Dye
    transition.

    (RSf, RSi, Ef, Ei, tf)
    """
    RSi = RS(2, 'F', 7 * 0.5)
    RSf = RS(2, 'G', 9 * 0.5)
    Ef = 154181.4942
    Ei = 170530.4040
    #Just the one transition
    #tf = 1.0 / 2e7
    #The mean:
    tf = 1.0/1.187e8
    return (RSf, RSi, Ef, Ei, tf)

def Diode_Values():
    """
    Returns the RS terms, energy levels, and lifetimes of the standard Diode
    transition.

    (RSf, RSi, Ef, Ei, tf)
    """
    RSi = RS(4, 'D', 5 * 0.5)
    RSf = RS(4, 'F', 7 * 0.5)
    Ei = 157673.4134
    Ef = 142717.0965
    #Just the one transition
    #tf = 1.0 / 1.07e7
    #The mean:
    tf = 1.0/1.326e8
    return (RSf, RSi, Ef, Ei, tf)

def ZeemanSpec_Padded(wld, s, RSi, RSf, Ei, Ef, ts, P, B, roll=0):
    """
        Creates a padded Zeeman array with 10x + 1 the number of elements of
        the original wavelength array, finds a spectrum, then returns the
        spectrum with the original number of points. This is so that
        undersampling the Zeeman spectrum is not a problem in deconvolution.
        INPUTS:
        wld     -   Wavelengths. The first and last values are used to generate
            the padded array, while the length is used to cut down the spectrum.
        s       -   Spectrum. Used to scale the zeeman profile.
        RSi     -   RS term for upper level of transition.
        RSf     -   RS term for lower level of transition.
        Ei      -   Energy for upper level of transition.
        Ef      -   Energy for lower level of transition.
        ts      -   Tuple of mean lifetimes for upper and lower states. A one
                    element tuple is fine.
        P       -   Polarization of desired transition.
        roll    -   (Optional) # to roll the padded wavelength array.
        OUTPUTS:
        (ZPSC, wldp, ZP)
        ZPSC    -   Zeeman spectrum reduced back to wld size.
        wldp    -   Padded wavelength array.
        ZP      -   Padded Zeeman profile.
    """
    wldp = np.linspace(wld[0], wld[-1], 10 * np.size(wld) + 1)
    lns = Zeeman_Lines(RSi, RSf, Ei, Ef, wldp, P, B, ts)
    lns = lns / np.max(lns) * np.max(s)
    [f, gl] = spec.spec(lns, wldp[1] - wldp[0])
    gl = gl * np.sqrt(float(np.size(gl)) / float(np.size(wld)))
    glc = gl[(np.size(gl)) / 2 - (np.size(wld)) / 2 : (np.size(gl))
            / 2 + (np.size(wld)) / 2 + 1]
    return (glc, wldp, lns)



def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
