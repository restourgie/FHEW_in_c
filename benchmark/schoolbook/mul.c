#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mul.h"
#include "fft/normal_fft.h"
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
* TANGENT FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void split_radix_mul(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[ZEROPAD];
  double complex cplx_y[ZEROPAD];
  double complex cplx_res[ZEROPAD];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
    cplx_x[i+REALDIM] = 0.0;
    cplx_y[i+REALDIM] = 0.0;
  }

  tangent_forward(cplx_x,ZEROPAD,0);

  // tangent_4_forward(cplx_y,ZEROPAD,0);

  // for (int i = 0; i < ZEROPAD; ++i)
  // {
  //   cplx_res[i] = (cplx_x[i] * cplx_y[i])/ZEROPAD;
  // }
  // Tangent_inverse(cplx_res,ZEROPAD,0);

  // for (int i = 0; i < REALDIM; ++i)
  // {
  //   r->v[i] = cplx_res[i] - cplx_res[i+REALDIM];
  // }
}

/******************************************************************
*
* SPLIT RADIX FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void split_radix_mul(ring_t *r, const ring_t *x, const ring_t *y)
{	
  double complex cplx_x[ZEROPAD];
  double complex cplx_y[ZEROPAD];
  double complex cplx_res[ZEROPAD];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
    cplx_x[i+REALDIM] = 0.0;
    cplx_y[i+REALDIM] = 0.0;
  }

  split_radix_recursive(cplx_x,ZEROPAD,0);
  split_radix_recursive(cplx_y,ZEROPAD,0);

  for (int i = 0; i < ZEROPAD; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/ZEROPAD;
  }
  split_radix_recursive_inverse(cplx_res,ZEROPAD,0);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = cplx_res[i] - cplx_res[i+REALDIM];
  }
}

/******************************************************************
*
* Twisted FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void twisted_mul(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[ZEROPAD];
  double complex cplx_y[ZEROPAD];
  double complex cplx_res[ZEROPAD];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
    cplx_x[i+REALDIM] = 0.0;
    cplx_y[i+REALDIM] = 0.0;
  }

  twisted_recursive(cplx_x,ZEROPAD,0);
  twisted_recursive(cplx_y,ZEROPAD,0);

  for (int i = 0; i < ZEROPAD; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/ZEROPAD;
  }
  twisted_recursive_inverse(cplx_res,ZEROPAD,0);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = cplx_res[i] - cplx_res[i+REALDIM];
  }
}

/******************************************************************
*
* SCHOOLBOOK NEGACYCLIC FFT
*
******************************************************************/
void normal_fft_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[ZEROPAD];
  double complex cplx_y[ZEROPAD];
  double complex cplx_res[ZEROPAD];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
    cplx_x[i+REALDIM] = 0.0;
    cplx_y[i+REALDIM] = 0.0;
  }

  recursive_FFT(cplx_x,ZEROPAD,0,1);
  recursive_FFT(cplx_y,ZEROPAD,0,1);

  for (int i = 0; i < ZEROPAD; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/ZEROPAD;
  }
  inverse_FFT(cplx_res,ZEROPAD,0,1);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = cplx_res[i] - cplx_res[i+REALDIM];
  }
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

