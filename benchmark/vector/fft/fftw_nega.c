#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include "../mul.h"
#include "fftw_nega.h"

fftw_complex *in,*out;
fftw_plan plan_fft_back,plan_fft_forw;

double complex table[CPLXDIM];



void int_table(){
  // table = malloc(CPLXDIM * sizeof(double complex));
  for (int i = 0; i < CPLXDIM; ++i)
  {
    table[i] = W(ROOTDIM,i);
  }
}

void FFTW_nega_setup() {
  in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * CPLXDIM);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * CPLXDIM);
  plan_fft_forw = fftw_plan_dft_1d(CPLXDIM, in, out, FFTW_FORWARD, FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_1d(CPLXDIM, out, in, FFTW_BACKWARD, FFTW_PATIENT);
  int_table();
}

/******************************************************************
*
* LOOKUPTABLE TWIST
*
******************************************************************/
void table_twist()
{
  for (int i = 1; i < CPLXDIM; ++i)
    in[i] = in[i] * table[i];
}

void table_untwist()
{
  for (int i = 1; i < CPLXDIM; ++i)
    in[i] = in[i] * conj(table[i]);

}

  
void FFTW_nega_forward(double complex *res, const ring_t *val) {
  for (int k = 0; k < CPLXDIM; ++k)
    in[k] = val->v[k] + I* val->v[CPLXDIM+k]; 
  table_twist(); 
  fftw_execute(plan_fft_forw); 
  for (int k = 0; k < CPLXDIM; ++k) 
    res[k] = (double complex) out[k];				
}

void FFTW_nega_backward(ring_t *res, const double complex *val){
  for (int k = 0; k < CPLXDIM; ++k) {
    out[k] = val[k];
  }
  fftw_execute(plan_fft_back);
  table_untwist();
  for (int k = 0; k < CPLXDIM; ++k){	
    res->v[k] = (long int) round(creal(in[k]));
    res->v[k+CPLXDIM] = (long int) round(cimag(in[k]));
  }
}