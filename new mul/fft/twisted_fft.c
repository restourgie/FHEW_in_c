#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "twisted_fft.h"
#include "support.h"

/******************************************************************
*
* TWISTED FFT MULTIPLICATION
*
******************************************************************/
void twisted_recursive_FFT(double complex *x,int n,int lo)
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
    twisted_recursive_FFT(x,m,lo);
    // printf("\n\n**************TWISTING**************\n");
    // print_complex(x,REALDIM);
    twist(x,n,m,lo+m);
    // print_complex(x,REALDIM);
    twisted_recursive_FFT(x,m,lo+m);
  }
}

void twisted_recursive_inverse_FFT(double complex *x,int n,int lo)
{
  if(n > 1)
  {
    int m = n/2;
    double complex temp;

    twisted_inverse_FFT(x,m,lo);
    twisted_inverse_FFT(x,m,lo+m);
    untwist(x,n,m,lo+m);

    for (int i = lo; i < lo+m; ++i)
    {
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }

  }
}