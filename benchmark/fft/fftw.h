#ifndef FFT_H
#define FFT_H

void FFTsetup();
void FFTWforward(double complex *res, const ring_t *val);
void FFTWbackward(ring_t *res, const double complex *val);
void FFTW_nega_forward(double complex *res, const ring_t *val);
void FFTW_nega_backward(ring_t *res, const double complex *val);

#endif