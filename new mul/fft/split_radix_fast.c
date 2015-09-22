// #include <xmmintrin.h>
#include <stdio.h>
#include <math.h>
#include "../support.h"

#define calc_cos(N,k) (cos(2.0 * M_PI * (double)k / (double) N))
#define calc_sin(N,k) (sin(2.0 * M_PI * (double)k / (double) N))

void fast_twist(cplx *cplx_x,int n,int m,int lo)
{
  double a,b,c,d;
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {	
  	a = cplx_x->real[i];
  	b = cplx_x->imag[i];	
  	c = calc_cos(n,j);
  	d = calc_sin(n,j);
  	 
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
  	cplx_x->real[i] = (a*c) - (b*d);
  	cplx_x->imag[i] = (a*d) + (b*c);
    ++j;
  }
}

void fast_untwist(cplx *cplx_x,int n,int m,int lo)
{
  double a,b,c,d;
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    a = cplx_x->real[i];
  	b = cplx_x->imag[i];	
  	c = calc_cos(n,j);
  	d = -calc_sin(n,j);
  	 
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
  	cplx_x->real[i] = (a*c) - (b*d);
  	cplx_x->imag[i] = (a*d) + (b*c);
    ++j;
  }
}

void split_radix_fast(cplx *x,int n,int lo)
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
		//Go from (x^4n +1 to x^2n -1 and x^2n +1)
		for(int i=lo; i < lo+m;++i){
		  temp_real = x->real[i];
		  temp_im = x->imag[i];

		  x->real[i] = temp_real + x->real[i+m];
		  x->imag[i] = temp_im + x->imag[i+m];

		  x->real[i+m] = temp_real - x->real[i+m];
		  x->imag[i+m] = temp_im - x->imag[i+m];
		}
		//Do recursive step for (x^2n -1)
		split_radix_fast(x,m,lo);

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
		fast_twist(x,n,m,lo);
		split_radix_fast(x,m,lo);
		fast_untwist(x,n,m,lo+m);
		split_radix_fast(x,m,lo+m);
	}
}

void split_radix_fast_inverse(cplx *x,int n,int lo)
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
    split_radix_fast_inverse(x,m,lo+m);
    fast_twist(x,n,m,lo+m);
    split_radix_fast_inverse(x,m,lo);
    fast_untwist(x,n,m,lo);
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
    split_radix_fast_inverse(x,m,lo);
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