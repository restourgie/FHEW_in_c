#ifndef FFT_H
#define FFT_H

void FFTsetup();
void FFTWforward(double complex *res, const ring_t *val);
void FFTWbackward(ring_t *res, const double complex *val);
void FFTWforward_original(double complex *res, const ring_t *val);
void FFTWbackward_original(ring_t *res, const double complex *val);

#endif