#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "mul.h"
#include "support.h"
#include "fft/fft_negacyc.h"
#include "fft/split_radix_fft.h"
#include "fft/split_radix_fast.h"

/******************************************************************
*
* SPLIT RADIX FFT FAST MULTIPLICATION
*
******************************************************************/
void split_radix_fast_mul(ring_t *r, const ring_t *x, const ring_t *y){
	printf("\n\n**************split-radix FAST**************\n");
	double cplx_x[REALDIM];
	double cplx_y[REALDIM];
	double cplx_res[REALDIM];

	for (int i = 0; i < REALDIM; ++i)
	{
		cplx_x[i] = x->v[i];
		cplx_y[i] = y->v[i];
	}
	printf("\n\n**************Normal X**************\n");
	print_double(cplx_x,CPLXDIM);
	fast_twist(cplx_x,2*REALDIM,CPLXDIM,0);
	printf("\n\n**************X AFTER TWIST**************\n");
	print_double(cplx_x,CPLXDIM);

}

/******************************************************************
*
* SPLIT RADIX FFT NEGACYCLIC MULTIPLICATION
*
******************************************************************/
void split_radix_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);
  printf("\n\n**************Normal X**************\n");
  print_complex(cplx_x,CPLXDIM);
  twist(cplx_x,2*REALDIM,CPLXDIM,0);
  printf("\n\n**************TWISTED X**************\n");
  print_complex(cplx_x,CPLXDIM);
  twist(cplx_y,2*REALDIM,CPLXDIM,0);

  split_radix_recursive(cplx_x,CPLXDIM,0);
  // printf("\n\n**************FFT X**************\n");
  // print_complex(cplx_x,CPLXDIM);

  split_radix_recursive(cplx_y,CPLXDIM,0);

  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }

  split_radix_recursive_inverse(cplx_res,CPLXDIM,0);
  untwist(cplx_res,2*REALDIM,CPLXDIM,0);
  // printf("\n\n**************split-radix FFT MUL RESULT**************\n");
  // print_complex(cplx_res,CPLXDIM);

  to_real(cplx_res,r);
}

void normal_fft_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
	smart_complex_mul(r,x,y);
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
    // printf("%f + i%f\n",creal(big[i+512]),cimag(big[i+512]));
    t = big[i] * (0+I*1);
    // printf("%f + i%f\n",creal(t),cimag(t));
    cplx_res[i-CPLXDIM] = big[i-CPLXDIM] + t;
  }
  // printf("\n\n**************NAIVE COMPLEX CALC**************\n");
  // print_complex(cplx_res,CPLXDIM);
  to_real(cplx_res,r);   
}
