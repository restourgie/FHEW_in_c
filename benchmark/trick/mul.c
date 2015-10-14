#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mul.h"
#include "fft/fft_negacyc.h"
#include "fft/split_radix_fft.h"
#include "fft/twisted_fft.h"
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
* SPLIT RADIX FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void split_radix_mul(ring_t *r, const ring_t *x, const ring_t *y)
{	
  // printf("*********************STARTING SPLIT RADIX***********************\n");
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);

  fft_sr_forward(cplx_x);
  // print_complex(cplx_x,CPLXDIM);
  fft_sr_forward(cplx_y);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  fft_sr_backward(cplx_res);
  // print_complex(cplx_res,CPLXDIM);
  to_real(cplx_res,r);
}

/******************************************************************
*
* Twisted FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void twisted_mul(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);

  fft_twisted_forward(cplx_x);
  fft_twisted_forward(cplx_y);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  fft_twisted_backward(cplx_res);
  to_real(cplx_res,r);
}

/******************************************************************
*
* SCHOOLBOOK NEGACYCLIC FFT
*
******************************************************************/
void normal_fft_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);

  double complex root = I;
  root = csqrt(root);
  // print_complex(cplx_x,CPLXDIM);
  recursive_phi(cplx_x,CPLXDIM,0,root);
  // print_complex(cplx_x,CPLXDIM);
  recursive_phi(cplx_y,CPLXDIM,0,root);


  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }

  inverse_phi(cplx_res,CPLXDIM,0,root);

  to_real(cplx_res,r);
}

/******************************************************************
*
* COMPLEX MULTIPLICATION
*
******************************************************************/
void naive_complex_mul(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);

  double complex t;
  double complex big[REALDIM];

  for(int i=0;i<REALDIM;++i)
    big[i] = 0;

  for(int i =0;i<CPLXDIM;++i)
    for(int j=0;j<CPLXDIM;++j)
      big[i+j] += cplx_x[i] * cplx_y[j];    

  for(int i=CPLXDIM;i<REALDIM;++i){
    t = big[i] * (0+I*1);
    cplx_res[i-CPLXDIM] = big[i-CPLXDIM] + t;
  }

  to_real(cplx_res,r);   
}
