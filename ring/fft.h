#ifndef FFT_H
#define FFT_H

#include <stdint.h>
#include <complex.h>

void FFTforward(double complex *r, const uint32_t *x);
void FFTbackward(uint32_t *r,  const double complex *x);
void CPLX_FFTforward(double complex *r);
void CPLX_FFTbackward(double complex *r);
// void FFTbackward(uint32_t *r, double complex *x);
void FFTsetup();

#endif
