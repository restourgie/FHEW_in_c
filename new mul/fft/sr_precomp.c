#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "../mul.h"

/******************************************************************
*
* SPLIT RADIX FFT MULTIPLICATION
*
******************************************************************/
void sr_precomp(cplx *x,int n,int lo)
{
  // print_double(x,CPLXDIM);
  double temp_real,temp_im;
  if(n == 2){
    temp_real = x->real[lo];
    temp_im = x->imag[lo];

    x->real[lo] = temp_real + x->real[lo+1];
    x->imag[lo] = temp_im + x->imag[lo+1];

    x->real[lo+1] = temp_real - x->real[lo+1];
    x->imag[lo+1] = temp_im - x->imag[lo+1];
  }
  else if(n > 2){
    int m = n/2;
    //Go from (x^4n-1 to x^2n -1 and x^2n +1)
    for(int i=lo; i < lo+m;++i){
      temp_real = x->real[i];
      temp_im = x->imag[i];

      x->real[i] = temp_real + x->real[i+m];
      x->imag[i] = temp_im + x->imag[i+m];

      x->real[i+m] = temp_real - x->real[i+m];
      x->imag[i+m] = temp_im - x->imag[i+m];
    }
    //Do recursive step for (x^2n -1)
    sr_precomp(x,m,lo);

    lo = lo+m;
    m = m/2;

    //Go from (x^2n +1 to x^n -i and x^n +i)
    for (int i = lo; i < lo+m; ++i)
    {
      temp_real = x->real[i];
      temp_im = x->imag[i];

      //(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
      x->real[i] = temp_real - x->imag[i+m];
      x->imag[i] = temp_im + x->real[i+m];

      //(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
      temp_real = temp_real + x->imag[i+m];
      x->imag[i+m] = temp_im - x->real[i+m];
      x->real[i+m] = temp_real;
    }
    table_twist(x,n,m,lo);
    sr_precomp(x,m,lo);
    table_untwist(x,n,m,lo+m);
    sr_precomp(x,m,lo+m);
  }
}

void sr_precomp_inverse(cplx *x,int n,int lo)
{
  double temp_real,temp_im;
  if(n == 2){
    temp_real = x->real[lo];
    temp_im = x->imag[lo];

    x->real[lo] = temp_real + x->real[lo+1];
    x->imag[lo] = temp_im + x->imag[lo+1];

    x->real[lo+1] = temp_real - x->real[lo+1];
    x->imag[lo+1] = temp_im - x->imag[lo+1];
  }
  else if(n > 2){
    // printf("n = %d lo = %d\n",n,lo );
    int m = n/4;
    lo = lo+n/2;
    // printf("m = %d lo = %d\n",m,lo );
    sr_precomp_inverse(x,m,lo+m);
    table_twist(x,n,m,lo+m);
    sr_precomp_inverse(x,m,lo);
    table_untwist(x,n,m,lo);
    // printf("BEFORE\n");
    // print_double(x,CPLXDIM);
    for (int i = lo; i < lo+m; ++i)
    {
    // printf("i = %d, m = %d\n",i,m );
      temp_real = x->real[i];
      temp_im = x->imag[i];

      x->real[i] = temp_real + x->real[i+m];
      x->imag[i] = temp_im + x->imag[i+m];

      temp_real = temp_real - x->real[i+m];
      temp_im = temp_im - x->imag[i+m];

    x->real[i+m] = temp_im;
    x->imag[i+m] = -temp_real;
    }
    // printf("AFTER\n");
    // print_double(x,CPLXDIM);
    m = m*2;
    lo = lo -m;
    // printf("m = %d lo = %d\n",m,lo );
    sr_precomp_inverse(x,m,lo);
    for(int i=lo; i < lo+m;++i){
      temp_real = x->real[i];
      temp_im = x->imag[i];

      x->real[i] = temp_real + x->real[i+m];
      x->imag[i] = temp_im + x->imag[i+m];

      x->real[i+m] = temp_real - x->real[i+m];
      x->imag[i+m] = temp_im - x->imag[i+m];
    }
  }
}