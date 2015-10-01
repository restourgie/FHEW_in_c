#ifndef FFT_NEGACYC_H
#define FFT_NEGACYC_H

#include "../mul.h"

void inverse_phi(double complex *x,int n,int lo,double complex root);
void recursive_phi(double complex *x,int n,int lo,double complex root);

#endif