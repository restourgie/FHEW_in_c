#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include "../mul.h"
#include "fftw.h"

double *in;
double complex table[CPLXDIM];
fftw_complex *out,*nega_in,*nega_out;
fftw_plan plan_fft_forw, plan_fft_back,plan_fft_nega_forw,plan_fft_nega_back;

void int_table(){
  // table = malloc(CPLXDIM * sizeof(double complex));
  for (int i = 0; i < CPLXDIM; ++i)
  {
    table[i] = W(ROOTDIM,i);
  }
} 

void FFTsetup() {
  in = (double*) fftw_malloc(sizeof(double) * 2*REALDIM);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (REALDIM + 2));
  plan_fft_forw = fftw_plan_dft_r2c_1d(2*REALDIM, in, out,  FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_c2r_1d(2*REALDIM, out, in,  FFTW_PATIENT);

  nega_in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * CPLXDIM);
  nega_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * CPLXDIM);
  plan_fft_nega_forw = fftw_plan_dft_1d(CPLXDIM, nega_in, nega_out, FFTW_FORWARD, FFTW_PATIENT);
  plan_fft_nega_back = fftw_plan_dft_1d(CPLXDIM, nega_out, nega_in, FFTW_BACKWARD, FFTW_PATIENT);
  int_table();
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

/******************************************************************
*
* LOOKUPTABLE TWIST
*
******************************************************************/
void table_twist()
{
  for (int i = 1; i < CPLXDIM; ++i)
    nega_in[i] = nega_in[i] * table[i];
}

void table_untwist()
{
  for (int i = 1; i < CPLXDIM; ++i)
    nega_in[i] = nega_in[i] * conj(table[i]);

}

  
void FFTW_nega_forward(double complex *res, const ring_t *val) {
  for (int k = 0; k < CPLXDIM; ++k)
    nega_in[k] = val->v[k] + I* val->v[CPLXDIM+k]; 
  table_twist(); 
  fftw_execute(plan_fft_nega_forw); 
  for (int k = 0; k < CPLXDIM; ++k) 
    res[k] = (double complex) nega_out[k];       
}

void FFTW_nega_backward(ring_t *res, const double complex *val){
  for (int k = 0; k < CPLXDIM; ++k) {
    nega_out[k] = val[k];
  }
  fftw_execute(plan_fft_nega_back);
  table_untwist();
  for (int k = 0; k < CPLXDIM; ++k){  
    res->v[k] = (long int) round(creal(in[k]));
    res->v[k+CPLXDIM] = (long int) round(cimag(in[k]));
  }
}