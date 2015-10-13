#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "twisted_fft.h"
#include "../mul.h"

void twisted_twist(double complex *cplx_x,int n,int m,int lo)
{
  // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    // printf("i = %d, j = %d\n",i,j);
    cplx_x[i] = cplx_x[i] * W(n,j);
    ++j;
  }
}

void twisted_untwist(double complex *cplx_x,int n,int m,int lo)
{
  // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    // printf("i = %d, j = %d\n",i,j);
    cplx_x[i] = cplx_x[i] * conj(W(n,j));
    ++j;
  }
}

/******************************************************************
*
* TWISTED FFT MULTIPLICATION
*
******************************************************************/
void twisted_recursive(double complex *x,int n,int lo)
{
  if(n > 1)
  {
    int m = n/2;
    double complex temp;

    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    // printf("n = %d\n",n );
    // printf("\n\n**************RECURSION**************\n");
    twisted_recursive(x,m,lo);
    // printf("\n\n**************TWISTING**************\n");
    // print_complex(x,REALDIM);
    twisted_twist(x,n,m,lo+m);
    // print_complex(x,REALDIM);
    twisted_recursive(x,m,lo+m);
  }
}

void twisted_recursive_inverse(double complex *x,int n,int lo)
{
  if(n > 1)
  {
    int m = n/2;
    double complex temp;

    twisted_recursive_inverse(x,m,lo);
    twisted_recursive_inverse(x,m,lo+m);
    twisted_untwist(x,n,m,lo+m);

    for (int i = lo; i < lo+m; ++i)
    {
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }

  }
}

void fft_twisted_forward(double complex *x)
{
  twisted_twist(x,ROOTDIM,CPLXDIM,0);
  twisted_recursive(x,CPLXDIM,0);
}

void fft_twisted_backward(double complex *x)
{
  twisted_recursive_inverse(x,CPLXDIM,0);
  twisted_untwist(x,ROOTDIM,CPLXDIM,0);
}