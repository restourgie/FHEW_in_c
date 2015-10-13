#ifndef FFT_NEGACYC_LUT_H
#define FFT_NEGACYC_LUT_H

#include <complex.h>

void init_negacyc();
void inverse_phi_lut(double complex *x,int n,int lo,int level);
void recursive_phi_lut(double complex *x,int n,int lo,int level);

#endif