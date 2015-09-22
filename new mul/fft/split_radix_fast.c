// #include <xmmintrin.h>
#include <stdio.h>
#include <math.h>
#include "../support.h"

#define calc_cos(N,k) (cos(2.0 * M_PI * (double)k / (double) N))
#define calc_sin(N,k) (sin(2.0 * M_PI * (double)k / (double) N))

void fast_twist(double *cplx_x,int n,int m,int lo)
{
  double a,b,c,d;
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {	
  	a = cplx_x[i];
  	b = cplx_x[i+CPLXDIM];	
  	c = calc_cos(n,j);
  	d = calc_sin(n,j);
  	 
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
  	cplx_x[i] = (a*c) - (b*d);
  	cplx_x[i+CPLXDIM] = (a*d) + (b*c);
    ++j;
  }
}

void split_radix_fast(double complex *x,int n,int lo)
{
  double temp_real,temp_im;
  if(n == 2){
    temp_real = x[lo];
    temp_im = x[lo+CPLXDIM];
    
    x[lo] = temp_real + x[lo+1];
    x[lo+CPLXDIM] = temp_im + x[lo+1+CPLXDIM];
    
    x[lo+1] = temp_real - x[lo+1];
    x[lo+1+CPLXDIM] = temp_im - x[lo+1+CPLXDIM];
  }
  else if(n > 2){
    int m = n/2;
    //Go from (x^4n +1 to x^2n -1 and x^2n +1)
    for(int i=lo; i < lo+m;++i){
      temp_real = x[i];
      temp_im = x[i+CPLXDIM];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    //Do recursive step for (x^2n -1)
    split_radix_recursive(x,m,lo);

    lo = lo+m;
    m = m/2;

    //Go from (x^2n +1 to x^n -i and x^n +i)
    for (int i = lo; i < lo+m; ++i)
    {
      temp = x[i];
      x[i] = temp + I * x[i+m];
      x[i+m] = temp - I * x[i+m];
    }
    twist(x,n,m,lo);
    split_radix_recursive(x,m,lo);
    untwist(x,n,m,lo+m);
    split_radix_recursive(x,m,lo+m);
  }
}

void split_radix_fast_inverse(double complex *x,int n,int lo)
{

}