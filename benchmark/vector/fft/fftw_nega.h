#ifndef FFT_NEGA_H
#define FFT_NEGA_H

void FFTW_nega_setup();
void FFTW_nega_forward(double complex *res, const ring_t *val);
void FFTW_nega_backward(ring_t *res, const double complex *val);

#endif