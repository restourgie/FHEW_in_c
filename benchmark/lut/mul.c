#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mul.h"
#include "fft/sr_precomp.h"
#include "fft/fft_negacyc_lut.h"
#include "fft/tangent_fft.h"

/******************************************************************
*
* SUPPORT CODE
*
******************************************************************/
void print_complex(const double complex *a, int N){
    for(int i=0;i<N;++i)
      printf("cplxpoly[%d] = %f + i * %f\n",i,creal(a[i]),cimag(a[i]));
    printf("\n");
}

void print_double(const cplx *x,int N){
  for (int i = 0; i < N; ++i)
    printf("cplxpoly[%d] = %f + i * %f\n",i,x->real[i],x->imag[i]);
  printf("\n");
}

/******************************************************************
*
* CONVERSION
*
******************************************************************/
void to_complex(const ring_t *x, double complex *cplx_x)
{
  for(int i=0;i<CPLXDIM;++i)
    cplx_x[i] = x->v[i] + I*x->v[i+CPLXDIM];
}

void to_real(const double complex *cplx_x, ring_t *x)
{
  for(int i=0;i<CPLXDIM;++i){
    x->v[i] = round(creal(cplx_x[i]));
    x->v[i+CPLXDIM] = round(cimag(cplx_x[i]));
  }
}

void init(){
  //PRECOMP TABLES FOR PRECOMP FFT
  init_table();
  //INIT LOOKUPTABLES NEGACYCLIC FFT
  init_negacyc();
  //INIT LOOKUPTABLES TANGENT FFT
  init_tangent();
}

/******************************************************************
*
* TANGENT FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void tangent_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
// { printf("*********************STARTING TANGENT***********************\n");
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);

  tangent_forward(cplx_x);
  // print_complex(cplx_x,CPLXDIM);
  tangent_forward(cplx_y);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  tangent_backward(cplx_res);
  // print_complex(cplx_res,CPLXDIM);
  to_real(cplx_res,r);
}

/******************************************************************
*
* SPLIT RADIX PRECOMPUTED FFT MULTIPLICATION
*
******************************************************************/
void sr_precomp_mul(ring_t *r, const ring_t *x, const ring_t *y){
  cplx cplx_x,cplx_y,cplx_res;

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x.real[i] = x->v[i];
    cplx_y.real[i] = y->v[i];
    cplx_x.imag[i] = x->v[j];
    cplx_y.imag[i] = y->v[j];
    ++j;
  }
  
  fft_precompsr_forward(&cplx_x);
  fft_precompsr_forward(&cplx_y);

  double a,b,c,d;
  for (int i = 0; i < CPLXDIM; ++i)
  {
      a = cplx_x.real[i];
      b = cplx_x.imag[i]; 
      c = cplx_y.real[i];
      d = cplx_y.imag[i];
      //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
      cplx_res.real[i] = ((a*c) - (b*d))/CPLXDIM;
      cplx_res.imag[i] = ((a*d) + (b*c))/CPLXDIM;
  }
  fft_precompsr_backward(&cplx_res,r);
}

/******************************************************************
*
* NEGACYCLIC FFT LOOK UP TABLE
*
******************************************************************/
void negacyc_lut_fft_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);
  // printf("*****************START LUT FFT********************\n");
  recursive_phi_lut(cplx_x,CPLXDIM,0,0);
  // print_complex(cplx_x,CPLXDIM);
  recursive_phi_lut(cplx_y,CPLXDIM,0,0);


  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }

  inverse_phi_lut(cplx_res,CPLXDIM,0,0);

  to_real(cplx_res,r);
}

/******************************************************************
*
*	NAIVE SCHOOLBOOK MULTIPLICATION
*
******************************************************************/
/* Very simple schoolbook multiplication. Works. */
void naive_real_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  int i,j;
  for(i=0;i<REALDIM;i++)
    r->v[i] = 0;

  for(i=0;i<REALDIM;i++)
  {
    for(j=0;j<REALDIM;j++)
    {
      if(i+j < REALDIM)
        r->v[i+j] += x->v[i] * y->v[j];
      else
        r->v[i+j-REALDIM] -= x->v[i] * y->v[j];
    }
  }
}