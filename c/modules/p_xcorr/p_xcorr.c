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
#include <omp.h>

#include "p_xcorr.h"

void corr3_parallel(double *f1, double *f2, int n_f1, int n_f2, int tau1, int tau2, 
        int unbiased, double *xcorr) {
    /*
     *This program finds the three point correlation for [+/- tau1, +/- tau2].
     *
     *void corr3_parallel(double *f1, double *f2, int n_f1, int n_f2, int tau1, int tau2, 
     *                      int unbiased, double *xcorr)
     *
     *
     * INPUTS:
     *
     *  double *f1 - First data array to be correlated.
     *  double *f2 - Second data array to be correlated.
     *  int n_f1    - Size of first data array.
     *  int n_f2    - Size of second data array.
     *  int tau1    - Max displacement in both directions for f1.
     *  int tau2    - Max displacement in both directions for f2.
     *  int unbiased - Boolean choosing unbiased (unbiased = 1) or unnormalized
     *  xcorr.
     *
     *  OUTPUTS:
     *
     *  double *xcorr  - 2D array of 3 point correlations. xcorr[0,0] is
     *                      [-tau1, -tau2] and xcorr[tau1 -1, tau2 -1] 
     *
     * The memory must be allocated for double *xcorr before this function is
     * run.
     *
     * It is recommended to reshape the output xcorr according to (2 * tau1 +
     * 1, 2 * tau2 + 1).
     *
     * This current iteration assumes that n_f1 = n_f2. This may be changed in
     * the future.
     *
     */
    int i=0, j=0, nthreads, l1, l2;
    printf("tau1 is %d, tau2 is %d\n", tau1, tau2);
    if (unbiased) {
        #pragma omp parallel
        {
            nthreads = omp_get_num_threads();
            printf("Number of threads is %d\n",nthreads);
            l1 = (tau1 >= n_f1) ? n_f1 - 1: tau1;
            l2 = (tau2 >= n_f2) ? n_f2 - 1: tau2;
            #pragma omp for private(i, j)
            for (i=0; i <= 2 * l1; i++) {
                for (j=0; j<= 2 * l2; j++) {
                    xcorr[i + j * (2 * l1 + 1)] = xcorr_unbiased(f1, f2, n_f1, i - l1, j - l2);
                }
            }
        }
    }else  {
        #pragma omp parallel
        {
            l1 = (tau1 >= n_f1) ? n_f1 - 1: tau1;
            l2 = (tau2 >= n_f2) ? n_f2 - 1: tau2;
            #pragma omp for private(i, j)
            for (i=0; i <= 2 * l1; i++) {
                for (j=0; j<= 2 * l2; j++) {
                    xcorr[i + j * (2 * l1 + 1)] = xcorr_sum(f1, f2, n_f1, i - l1, j - l2);
                }
            }
        }
    }
    return;
}

double xcorr_sum(double *f1, double *f2, int n, int tau1, int tau2) {
    //Returns the single element three point correlation at the point given by
    //tau1, tau2.

    int i, nf, ni;
    double sum = 0;

    //    for (i=0; i < n; i++) printf("f1 = %f, f2 = %f\n", f1[i], f2[i]);
    //if tau1 > tau2, make ni = tau1, else ni = tau2
    if (tau1 > 0 && tau2 > 0) {
        ni = (tau1 > tau2) ? tau1 : tau2;
        nf = n;
    }
    else if (tau1 <= 0 && tau2 <= 0)  {
        ni = 0;
        nf = (tau1 < tau2) ? n - abs(tau1) : n - abs(tau2);
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
    //    printf("ni = %d, nf = %d\n", ni, nf);
    //    printf("\nFor tau1 = %d, tau2 = %d, we are summing over:\n", tau1, tau2);
    for (i=ni; i < nf; i++) {
        sum += f1[i] * f1[i - tau1] * f2[i - tau2]; 
    }
    //    printf("Sum done. \n\n");
    //    printf("sum = %f\n", sum);
    return sum;
}

double xcorr_unbiased(double *f1, double *f2, int n, int tau1, int tau2) {
    //Returns the single element three point correlation at the point given by
    //tau1, tau2.
    int i, nf, ni;
    double sum = 0;

    //if tau1 > tau2, make ni = tau1, else ni = tau2
    if (tau1 > 0 && tau2 > 0) {
        ni = (tau1 > tau2) ? tau1 : tau2;
        nf = n;
    }
    else if (tau1 <= 0 && tau2 <= 0)  {
        ni = 0;
        nf = (tau1 < tau2) ? n - abs(tau1) : n - abs(tau2);
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
    sum = sum / (nf - ni);
    return sum;
}



