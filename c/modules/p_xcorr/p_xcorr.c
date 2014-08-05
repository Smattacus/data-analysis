/****
 *
 * This is a couple of standalone functions for calculating the 2nd order cross
 * correlation using the OpenMP (shared memory) parallel processing framework.
 *
 * The formula followed for the cross correlation is 
 *
 *  Corr3(tau1, tau2) = sum over t of (f_1(t) * f_1(t - tau1) * f_2(t - tau1)
 *
 * Where t is the index of the arrays for f_1(t) and f_(t).
 *
 * Sean Mattingly
 * University of Iowa
 * August 2014
 *
 * ***/
#include <stdio.h>
#include <stdlib.h>
#include "p_xcorr.h"

void corr3_parallel(double *f1, double *f2, int n_f1, int n_f2, int tau1, int tau2, double *xcorr) {
    /*
     *This program finds the three point correlation for [+/- tau1, +/- tau2].
     *INPUTS:
     *
     *  double *f1 - First data array to be correlated.
     *  double *f2 - Second data array to be correlated.
     *  int n_f1    - Size of first data array.
     *  int n_f2    - Size of second data array.
     *  int tau1    - Max displacement in both directions for f1.
     *  int tau2    - Max displacement in both directions for f2.
     *
     *  OUTPUTS:
     *
     *  double *xcorr  - 2D array of 3 point correlations. xcorr[0,0] is
     *                      [-tau1, -tau2] and xcorr[tau1 -1, tau2 -1] 
     *
     * The memory must be allocated for double *xcorr before this function is
     * run.
     *
     * It is recommended to reshape the output xcorr according to (2 * tau1 -
     * 1, 2 * tau2 - 1).
     *
     * This current iteration assumes that n_f1 = n_f2. This may be changed in
     * the future.
     *
     */
    int i=0, j=0;
    double temp;

    for (i=(-tau1+1); i < tau1; i++) {
        for (j=(-tau2 + 1); j < tau2; j++) {
            xcorr[i + tau1 + j * (2 * tau1 - 1)] = xcorr_sum(f1, f2, n_f1, i, j); 
            printf("Output is %d", xcorr[i + tau1 + j * (2 * tau1 -1));
        }
    }
    return;
}

double xcorr_sum(double *f1, double *f2, int n, int tau1, int tau2) {
    //Returns the single element three point correlation at the point given by
    //tau1, tau2.
    
    int i, nf, ni;
    double sum = 0;

    //if tau1 > tau2, make ni = tau1, else ni = tau2
    if (tau1 > 0 && tau2 > 0) ni = (tau1 > tau2) ? tau1 : tau2;
    else if (tau1 <= 0 && tau2 <= 0)  {
        ni = 0;
        nf = abs(tau1);
    }
    //One positive, one negative, set ni = positive tau
    else {
    //set ni = positive tau, nf = n - abs(negative tau)
        if (tau1 > tau2) {
            ni = tau1;
            nf = n - abs(tau2);
        } else {
            ni = tau2;
            nf = n - abs(tau1);
        }
    
    }

    for (i=ni; i < nf; i++) {
        sum += f1[i] * f1[i - tau1] * f2[i - tau2]; 
    }
    sum /= (n - ni + 1);
    return sum;
}


