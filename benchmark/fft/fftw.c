#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include "mul.h"
#include "fftw.h"

fftw_complex *nega_in,*nega_out;
fftw_plan plan_fft_nega_forw,plan_fft_nega_back;


void FFTsetup() {
  nega_in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * CPLXDIM);
  nega_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * CPLXDIM);
  plan_fft_nega_forw = fftw_plan_dft_1d(CPLXDIM, nega_in, nega_out, FFTW_FORWARD, FFTW_PATIENT);
  plan_fft_nega_back = fftw_plan_dft_1d(CPLXDIM, nega_out, nega_in, FFTW_BACKWARD, FFTW_PATIENT);
}

  
void FFTWforward(double complex *res, const ring_t *val) {
  for (int k = 0; k < CPLXDIM; ++k)
    nega_in[k] = val->v[k] + I* val->v[CPLXDIM+k];  
  fftw_execute(plan_fft_nega_forw); 
  for (int k = 0; k < CPLXDIM; ++k) 
    res[k] = (double complex) nega_out[k];       
}

void FFTWbackward(ring_t *res, const double complex *val){
  for (int k = 0; k < CPLXDIM; ++k) {
    nega_out[k] = val[k];
  }
  fftw_execute(plan_fft_nega_back);
  for (int k = 0; k < CPLXDIM; ++k){  
    res->v[k] = (long int) round(creal(nega_in[k]));
    res->v[k+CPLXDIM] = (long int) round(cimag(nega_in[k]));
  }
}