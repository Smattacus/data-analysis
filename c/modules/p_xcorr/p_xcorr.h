#include <omp.h>
#include <math.h>

double xcorr_sum(double *, double *, int, int, int);
void corr3_parallel(double *, double *, int, int, int, int, int, double *);
double xcorr_unbiased(double *, double *, int, int, int);
