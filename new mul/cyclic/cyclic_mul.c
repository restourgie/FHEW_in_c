#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "cyclic_mul.h"
#include "support.h"
#include "../fft/split_radix_fft.h"
#include "../fft/normal_fft.h"
#include "../fft/twisted_fft.h"

/******************************************************************
*
* NAIVE CYCLIC SCHOOLBOOK MULTIPLICATION
*
******************************************************************/
/* Very simple schoolbook multiplication. Works. */
void naive_cyclic_real_mul(ring_t *r, const ring_t *x, const ring_t *y)
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
        r->v[i+j-REALDIM] += x->v[i] * y->v[j];
    }
  }
}

/******************************************************************
*
* Normal FFT MULTIPLICATION
*
******************************************************************/
void normal_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[REALDIM];
  double complex cplx_y[REALDIM];
  double complex cplx_res[REALDIM];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
  }

  double complex root = 1;
  
  // printf("\n\n**************APPLY FFT**************\n");
  recursive_FFT(cplx_x,REALDIM,0,root);
  // print_complex(cplx_x,REALDIM);



  //printf("\n\n**************NORMAL FFT Y**************\n");
  recursive_FFT(cplx_y,REALDIM,0,root);
  //print_complex(cplx_y,CPLXDIM);


  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/REALDIM;
  }

  // printf("\n\n**************NORMAL FFT MUL RESULT**************\n");
  inverse_FFT(cplx_res,REALDIM,0,root);
  
  // print_complex(cplx_res,REALDIM);
  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = round(creal(cplx_res[i]));
  } 
}

/******************************************************************
*
* TWISTED FFT MULTIPLICATION
*
******************************************************************/
void twisted_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[REALDIM];
  double complex cplx_y[REALDIM];
  double complex cplx_res[REALDIM];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
  }
  // printf("\n\n**************APPLY TWISTED FFT**************\n");
  // print_complex(cplx_x,REALDIM);
  twisted_recursive_FFT(cplx_x,REALDIM,0);
  // print_complex(cplx_x,REALDIM);
  twisted_recursive_FFT(cplx_y,REALDIM,0);

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/REALDIM;
  }
  // printf("\n\n**************TWISTED FFT MUL RESULT Before inverse**************\n");
  // print_complex(cplx_res,REALDIM);
  printf("\n\n**************TWISTED FFT MUL RESULT**************\n");
  twisted_inverse_FFT(cplx_res,REALDIM,0);
  print_complex(cplx_res,REALDIM);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = round(creal(cplx_res[i]));
  }
}

/******************************************************************
*
* SPLIT RADIX FFT MULTIPLICATION
*
******************************************************************/
void split_radix_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[REALDIM];
  double complex cplx_y[REALDIM];
  double complex cplx_res[REALDIM];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
  }
  // printf("\n\n**************APPLY split-radix FFT**************\n");
  // print_complex(cplx_x,REALDIM);
  split_radix_recursive(cplx_x,REALDIM,0);

  // print_complex(cplx_x,REALDIM);
  split_radix_recursive(cplx_y,REALDIM,0);

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/REALDIM;
  }
  // printf("\n\n**************TWISTED FFT MUL RESULT Before inverse**************\n");
  // print_complex(cplx_res,REALDIM);
  // printf("\n\n**************split-radix FFT MUL RESULT**************\n");
  split_radix_recursive_inverse(cplx_res,REALDIM,0);
  // print_complex(cplx_res,REALDIM);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = round(creal(cplx_res[i]));
  }
}

