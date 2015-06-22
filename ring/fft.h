#ifndef FFT_H
#define FFT_H

#include <stdint.h>
#include <complex.h>

void FFTforward(double complex *r, const uint32_t *x);
void FFTbackward(uint32_t *r,  const double complex *x);

#endif
