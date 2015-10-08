#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include "fftw.h"
#include "../mul.h"

double *in;
fftw_complex *out;
fftw_plan plan_fft_forw, plan_fft_back;
  
void FFTsetup() {
  in = (double*) fftw_malloc(sizeof(double) * 2*REALDIM);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (REALDIM + 2));
  plan_fft_forw = fftw_plan_dft_r2c_1d(2*REALDIM, in, out,  FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_c2r_1d(2*REALDIM, out, in,  FFTW_PATIENT);
}
  
void FFTWforward(double complex *res, const ring_t *val) {
  for (int k = 0; k < REALDIM; ++k)	{
    in[k] = val->v[k];
    in[k+REALDIM] = 0.0;			
  }
  fftw_execute(plan_fft_forw); 
  for (int k = 0; k < CPLXDIM+1; ++k) 
    res[k] = (double complex) out[2*k+1];				
}

void FFTWbackward(ring_t *res, const double complex *val){
  for (int k = 0; k < CPLXDIM+1; ++k) {
    out[2*k+1] = (double complex) val[k]/REALDIM;
    out[2*k]   = (double complex) 0;
  }
  fftw_execute(plan_fft_back); 
  for (int k = 0; k < REALDIM; ++k)	
    res->v[k] = (long int) round(in[k]);
}