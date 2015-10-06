#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include "../mul.h"

void sr_vector(cplx_ptr *x,int n,int lo){
  // print_double(x,CPLXDIM);
  double temp_real,temp_im;
  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp;
  // __m128d real_small_x,imag_small_x,real_small_y,imag_small_y,imag_small_temp,real_small_temp;
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
    for(int i=lo; i < lo+m;i+=4){
      //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
      real_x = _mm256_load_pd(x->real+i);
      imag_x = _mm256_load_pd(x->imag+i);
      real_y = _mm256_load_pd(x->real+i+m);
      imag_y = _mm256_load_pd(x->imag+i+m);

      //TEMP IS X - Y 
      real_temp = _mm256_sub_pd(real_x,real_y);
      imag_temp = _mm256_sub_pd(imag_x,imag_y);

      //X = X+Y
      real_x = _mm256_add_pd(real_x,real_y);
      imag_x = _mm256_add_pd(imag_x,imag_y);

      _mm256_store_pd(x->real+i,real_x);
	  _mm256_store_pd(x->imag+i,imag_x);
	  _mm256_store_pd(x->real+i+m,real_temp);
	  _mm256_store_pd(x->imag+i+m,imag_temp);
    }
    //Do recursive step for (x^2n -1)
    sr_vector(x,m,lo);

    lo = lo+m;
    m = m/2;

    //Go from (x^2n +1 to x^n -i and x^n +i)
    for (int i = lo; i < lo+m; i+=4)
    {

      //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
      real_x = _mm256_load_pd(x->real+i);
      imag_x = _mm256_load_pd(x->imag+i);
      real_y = _mm256_load_pd(x->real+i+m);
      imag_y = _mm256_load_pd(x->imag+i+m);	  

      //(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
      real_temp = _mm256_add_pd(real_x,imag_y);
      imag_temp = _mm256_sub_pd(imag_x,real_y);

      //(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
      real_x = _mm256_sub_pd(real_x,imag_y);
      imag_x = _mm256_add_pd(imag_x,real_y);
      
      _mm256_store_pd(x->real+i,real_x);
	  _mm256_store_pd(x->imag+i,imag_x);
	  _mm256_store_pd(x->real+i+m,real_temp);
	  _mm256_store_pd(x->imag+i+m,imag_temp);
    }
    vector_twist(x,n,m,lo);
    sr_vector(x,m,lo);
    vector_untwist(x,n,m,lo+m);
    sr_vector(x,m,lo+m);
  }
}


void sr_vector_inverse(cplx_ptr *x,int n,int lo){
	
}