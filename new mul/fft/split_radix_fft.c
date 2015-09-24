#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "../mul.h"

// #define W(N,k) (cexp(2.0 * M_PI * I * (double)k / (double) N))

// void twist(double complex *cplx_x,int n,int m,int lo)
// {
//   // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
//   int j = 1;
//   for (int i = lo+1; i < lo+m; ++i)
//   {
//     // printf("i = %d, j = %d\n",i,j);
//     cplx_x[i] = cplx_x[i] * W(n,j);
//     ++j;
//   }
// }


// void untwist(double complex *cplx_x,int n,int m,int lo)
// {
//   // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
//   int j = 1;
//   for (int i = lo+1; i < lo+m; ++i)
//   {
//     // printf("i = %d, j = %d\n",i,j);
//     cplx_x[i] = cplx_x[i] * conj(W(n,j));
//     ++j;
//   }
// }

/******************************************************************
*
* SPLIT RADIX FFT MULTIPLICATION
*
******************************************************************/
void split_radix_recursive(double complex *x,int n,int lo)
{

  double complex temp;
  if(n == 2){
    temp = x[lo];
    x[lo] = temp + x[lo+1];
    x[lo+1] = temp - x[lo+1];
  }
  else if(n > 2){
    int m = n/2;
    //Go from (x^4n +1 to x^2n -1 and x^2n +1)
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    //Do recursive step for (x^2n -1)
    split_radix_recursive(x,m,lo);

    lo = lo+m;
    m = m/2;

    // printf("BEFORE THE I STEP\n");
    // print_complex(x,CPLXDIM);
    //Go from (x^2n +1 to x^n -i and x^n +i)
    for (int i = lo; i < lo+m; ++i)
    {
      // printf("i = %d, m = %d\n",i,m );
      temp = x[i];
      x[i] = temp + I * x[i+m];
      x[i+m] = temp - I * x[i+m];
    }
    // printf("AFTER THE I STEP\n");
    // print_complex(x,CPLXDIM);
    twist(x,n,m,lo);
    split_radix_recursive(x,m,lo);
    untwist(x,n,m,lo+m);
    // printf("AFTER UNTWISTING\n");
    // print_complex(x,CPLXDIM);
    split_radix_recursive(x,m,lo+m);
  }
}

void split_radix_recursive_inverse(double complex *x,int n,int lo)
{
  // printf("N = %d\n",n);
  // print_complex(x,REALDIM);
  double complex temp;
  if(n == 2){
    temp = x[lo];
    x[lo] = temp + x[lo+1];
    x[lo+1] = temp - x[lo+1];
  }
  else if(n > 2){
    // printf("n = %d lo = %d\n",n,lo );
    int m = n/4;
    lo = lo+n/2;
    // printf("m = %d lo = %d\n",m,lo );
    split_radix_recursive_inverse(x,m,lo+m);
    twist(x,n,m,lo+m);
    split_radix_recursive_inverse(x,m,lo);
    untwist(x,n,m,lo);
    // printf("BEFORE\n");
    // print_complex(x,CPLXDIM);
    for (int i = lo; i < lo+m; ++i)
    {
      // printf("i = %d, m = %d\n",i,m );
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = (temp - x[i+m])*-I;
    }
    // printf("AFTER\n");
    // print_complex(x,CPLXDIM);
    m = m*2;
    lo = lo -m;
    // printf("m = %d lo = %d\n",m,lo );
    split_radix_recursive_inverse(x,m,lo);
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
  }
}