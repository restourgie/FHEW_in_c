#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include "../mul.h"

void sr_vector(cplx_ptr *x,int n,int lo){
  double temp_real,temp_im;
  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp;
  __m128d real_x_small,imag_x_small,imag_y_small,real_y_small,imag_temp_small,real_temp_small;
  if(n > 8)
  {
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
    // printf("AFTER FIRST FOLD\n");
    //print_cplx(*x);
    //Do recursive step for (x^2n-1)
    sr_vector(x,m,lo);

    lo = lo+m;
    m = m/2;

    // printf("ENTERING m = %d\n",m);
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
    // printf("BEFORE VECTOR TWIST\n");
    // print_cplx(*x);
    vector_untwist(x,n,m,lo+m);
    // printf("AFTER VECTOR TWIST\n");
    // print_cplx(*x);
    sr_vector(x,m,lo+m);
  }
  else if(n == 8)
  {
    int m = 4;
    //Go from (x^8+1) to (x^4-1) and (x^4+1)
    //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
    real_x = _mm256_load_pd(x->real+lo);
    imag_x = _mm256_load_pd(x->imag+lo);
    real_y = _mm256_load_pd(x->real+lo+m);
    imag_y = _mm256_load_pd(x->imag+lo+m);

    //TEMP IS X - Y 
    real_temp = _mm256_sub_pd(real_x,real_y);
    imag_temp = _mm256_sub_pd(imag_x,imag_y);

    //X = X+Y
    real_x = _mm256_add_pd(real_x,real_y);
    imag_x = _mm256_add_pd(imag_x,imag_y);

    _mm256_store_pd(x->real+lo,real_x);
    _mm256_store_pd(x->imag+lo,imag_x);
    _mm256_store_pd(x->real+lo+m,real_temp);
    _mm256_store_pd(x->imag+lo+m,imag_temp);
    // printf("AFTER FIRST FOLD\n");
    //print_cplx(*x);
    //Do recursive step for (x^2n-1)
    sr_vector(x,m,lo);

    lo = lo+m;
    m = m/2;
    //Go from (x^4-1) to (x^2-i) and (x^2+i)
    //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
    real_x_small = _mm_load_pd(x->real+lo);
    imag_x_small = _mm_load_pd(x->imag+lo);
    real_y_small = _mm_load_pd(x->real+lo+m);
    imag_y_small = _mm_load_pd(x->imag+lo+m);   

    //(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
    real_temp_small = _mm_add_pd(real_x_small,imag_y_small);
    imag_temp_small = _mm_sub_pd(imag_x_small,real_y_small);

    //(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
    real_x_small = _mm_sub_pd(real_x_small,imag_y_small);
    imag_x_small = _mm_add_pd(imag_x_small,real_y_small);
    
    _mm_store_pd(x->real+lo,real_x_small);
    _mm_store_pd(x->imag+lo,imag_x_small);
    _mm_store_pd(x->real+lo+m,real_temp_small);
    _mm_store_pd(x->imag+lo+m,imag_temp_small);
    //WE ARE GOING TO DO A MANUAL TWIST SINCE 
    //(x->real[lo+1],x->imag[lo+1]) AND (x->real[lo+m+1],x->imag[lo+m+1])
    //ONLY NEED TO BE TWISTED
    int other_lo = lo+m;
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    temp_real = (x->real[lo+1] * LUT1[0][1]) - (x->imag[lo+1] * LUT2[0][1]);//ac -bd
    x->imag[lo+1] = (x->real[lo+1] * LUT2[0][1]) + (x->imag[lo+1] * LUT1[0][1]);
    x->real[lo+1] = temp_real;
    //OTHER TWIST
    temp_real = (x->real[other_lo+1] * LUT1[0][1]) - (x->imag[other_lo+1] * LUT3[0][1]);//ac -bd
    x->imag[other_lo+1] = (x->real[other_lo+1] * LUT3[0][1]) + (x->imag[other_lo+1] * LUT1[0][1]);
    x->real[other_lo+1] = temp_real;

    //WE CAN NOW FINISH THE FFT BY ADDING AND SUBTRACTING
    //LEFT PART
    temp_real = x->real[lo];
    temp_im = x->imag[lo];
    x->real[lo] = temp_real + x->real[lo+1];
    x->imag[lo] = temp_im + x->imag[lo+1];
    x->real[lo+1] = temp_real - x->real[lo+1];
    x->imag[lo+1] = temp_im - x->imag[lo+1];
    //RIGHT PART
    temp_real = x->real[other_lo];
    temp_im = x->imag[other_lo];
    x->real[other_lo] = temp_real + x->real[other_lo+1];
    x->imag[other_lo] = temp_im + x->imag[other_lo+1];
    x->real[other_lo+1] = temp_real - x->real[other_lo+1];
    x->imag[other_lo+1] = temp_im - x->imag[other_lo+1];
  }
  else if(n == 4)
  {
    // printf("ENTERING\n");
    int m = 2;
    //Go from (x^4-1) to (x^2-1) and (x^2+1)
    real_x_small = _mm_load_pd(x->real+lo);
    imag_x_small = _mm_load_pd(x->imag+lo);
    real_y_small = _mm_load_pd(x->real+lo+m);
    imag_y_small = _mm_load_pd(x->imag+lo+m);

    real_temp_small = _mm_sub_pd(real_x_small,real_y_small);
    imag_temp_small = _mm_sub_pd(imag_x_small,imag_y_small);
    real_x_small = _mm_add_pd(real_x_small,real_y_small);
    imag_x_small = _mm_add_pd(imag_x_small,imag_y_small);

    _mm_store_pd(x->real+lo,real_x_small);
    _mm_store_pd(x->imag+lo,imag_x_small);
    _mm_store_pd(x->real+lo+m,real_temp_small);
    _mm_store_pd(x->imag+lo+m,imag_temp_small);
    //DO LEFT BRANCH GO FROM (x^2-1) to (x-1) and (x+1)
    temp_real = x->real[lo];
    temp_im = x->imag[lo];
    x->real[lo] = temp_real + x->real[lo+1];
    x->imag[lo] = temp_im + x->imag[lo+1];
    x->real[lo+1] = temp_real - x->real[lo+1];
    x->imag[lo+1] = temp_im - x->imag[lo+1];
    //DO RIGHT BRANCH GO FROM (x^2+1) to (x-i) and (x+i) 
    lo = lo+m;
    temp_real = x->real[lo];
    temp_im = x->imag[lo];
    //(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
    x->real[lo] = temp_real - x->imag[lo+1];
    x->imag[lo] = temp_im + x->real[lo+1];
    //(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
    temp_real = temp_real + x->imag[lo+1];
    x->imag[lo+1] = temp_im - x->real[lo+1];
    x->real[lo+1] = temp_real;
  }
}

void sr_vector_inverse(cplx_ptr *x,int n,int lo){
  double temp_real,temp_im;
  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp;
  __m128d real_x_small,imag_x_small,imag_y_small,real_y_small,imag_temp_small,real_temp_small;
  if(n > 8)
  {
    //WHEN n > 8 we have the general problem, in order to get (x^(4n)-1) we need to collect
    //The left branch (x^(2n)-1) and the two right branches (x^n-1) and (x^n-1)
    //We need to twist the right branches to (x^n-i) and (x^n+i) in order to fold them back to (x^(2n)+1)
    //Afterwards we can fold (x^(2n)-1) and (x^(2n)+1) to get (x^4n-1)
    int m = n/4;
    lo = lo+n/2;
    // printf("m = %d lo = %d\n",m,lo);
    //Get right branch (x^n-1)
    sr_vector_inverse(x,m,lo+m);
    //Twist right branch to get (x^n+i)
    vector_twist(x,n,m,lo+m);
    //Get left branch (x^n-1)
    sr_vector_inverse(x,m,lo);
    //unTwist left branch to get (x^n-i)
    vector_untwist(x,n,m,lo);
    // printf("BEFORE STITCHING\n");
    // print_cplx(*x,CPLXDIM);
    //TIE (x^n-i) and (x^n+i) back to (x^(2n)+1)
    for (int i = lo; i < lo+m; i+=4)
    {
      //real_temp = a-d, imag_temp = b+c, real_y = a+d, imag_y = b-c
      //2a = (a-d)+(a+d), 2b = (b+c)+(b-c), 2c = (b+c)-(b-c), 2d = (a+d)-(a-d) 
      //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
      real_temp = _mm256_load_pd(x->real+i);
      imag_temp = _mm256_load_pd(x->imag+i);
      real_y = _mm256_load_pd(x->real+i+m);
      imag_y = _mm256_load_pd(x->imag+i+m);   

      real_x = _mm256_add_pd(real_temp,real_y);
      imag_x = _mm256_add_pd(imag_temp,imag_y);

      real_temp = _mm256_sub_pd(real_y,real_temp);
      real_y = _mm256_sub_pd(imag_temp,imag_y);
      
      _mm256_store_pd(x->real+i,real_x);
      _mm256_store_pd(x->imag+i,imag_x);
      _mm256_store_pd(x->real+i+m,real_y);
      _mm256_store_pd(x->imag+i+m,real_temp);
    }
    // printf("AFTER STITCHING\n");
    // print_cplx(*x,CPLXDIM);
    m = m*2;
    lo = lo -m;
    //GET (x^(2n)-1)
    sr_vector_inverse(x,m,lo);
    //TIE (x^(2n)-1) and (x^(2n)+1) to get (x^4n-1)
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
  }
  if(n == 8)
  {
    int m   = 2;
    lo  = lo+4;
    int other_lo = lo+m; 
    //WE CAN NOW START THE FFT BY TIEING TOGETHER (x-1) and (x+1) to (x^2-1)
    //LEFT PART
    temp_real = x->real[lo];
    temp_im = x->imag[lo];
    x->real[lo] = temp_real + x->real[lo+1];
    x->imag[lo] = temp_im + x->imag[lo+1];
    x->real[lo+1] = temp_real - x->real[lo+1];
    x->imag[lo+1] = temp_im - x->imag[lo+1];
    //RIGHT PART
    temp_real = x->real[other_lo];
    temp_im = x->imag[other_lo];
    x->real[other_lo] = temp_real + x->real[other_lo+1];
    x->imag[other_lo] = temp_im + x->imag[other_lo+1];
    x->real[other_lo+1] = temp_real - x->real[other_lo+1];
    x->imag[other_lo+1] = temp_im - x->imag[other_lo+1];
    //NOW WE NEED TO DO A MANUAL TWIST TO GO FROM (x^2-1) to (x^2-i) and (x^2-1) to (x^2+i)
    //WE DO AN INVERSE TWIST ON THE LEFT PART TO GO FROM (x^2-1) to (x^2-i)
    temp_real = (x->real[lo+1] * LUT1[0][1]) - (x->imag[lo+1] * LUT3[0][1]);//ac -bd
    x->imag[lo+1] = (x->real[lo+1] * LUT3[0][1]) + (x->imag[lo+1] * LUT1[0][1]);
    x->real[lo+1] = temp_real;
    //WE DO AN FORWARD TWIST ON THE RIGHT PART TO GO FROM (x^2-1) to (x^2+i)
    temp_real = (x->real[other_lo+1] * LUT1[0][1]) - (x->imag[other_lo+1] * LUT2[0][1]);//ac -bd
    x->imag[other_lo+1] = (x->real[other_lo+1] * LUT2[0][1]) + (x->imag[other_lo+1] * LUT1[0][1]);
    x->real[other_lo+1] = temp_real;
    //WE NOW NEED TO TIE (x^2-i) and (x^2+i) TOGETHER TO FORM (X^4+1)
    //real_temp = a-d, imag_temp = b+c, real_y = a+d, imag_y = b-c
    //2a = (a-d)+(a+d), 2b = (b+c)+(b-c), 2c = (b+c)-(b-c), 2d = (a+d)-(a-d) 
    //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
    real_temp_small = _mm_load_pd(x->real+lo);
    imag_temp_small = _mm_load_pd(x->imag+lo);
    real_y_small = _mm_load_pd(x->real+lo+m);
    imag_y_small = _mm_load_pd(x->imag+lo+m);   

    real_x_small = _mm_add_pd(real_temp_small,real_y_small);
    imag_x_small = _mm_add_pd(imag_temp_small,imag_y_small);

    real_temp_small = _mm_sub_pd(real_y_small,real_temp_small);
    real_y_small = _mm_sub_pd(imag_temp_small,imag_y_small);
    
    _mm_store_pd(x->real+lo,real_x_small);
    _mm_store_pd(x->imag+lo,imag_x_small);
    _mm_store_pd(x->real+lo+m,real_y_small);
    _mm_store_pd(x->imag+lo+m,real_temp_small);
    //NOW WE HAVE THE ENTIRE RIGHT BRANCH (x^4+1)
    //RECURSIVE CALL TO GET LEFT BRANCH (x^4-1)
    m = 4;
    lo = lo-m;
    sr_vector_inverse(x,m,lo);
    //NOW WE NEED TO TIE TOGETHER (x^4-1) and (x^4+1) to form (x^8-1)
    //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
    real_x = _mm256_load_pd(x->real+lo);
    imag_x = _mm256_load_pd(x->imag+lo);
    real_y = _mm256_load_pd(x->real+lo+m);
    imag_y = _mm256_load_pd(x->imag+lo+m);

    //TEMP IS X - Y 
    real_temp = _mm256_sub_pd(real_x,real_y);
    imag_temp = _mm256_sub_pd(imag_x,imag_y);

    //X = X+Y
    real_x = _mm256_add_pd(real_x,real_y);
    imag_x = _mm256_add_pd(imag_x,imag_y);

    _mm256_store_pd(x->real+lo,real_x);
    _mm256_store_pd(x->imag+lo,imag_x);
    _mm256_store_pd(x->real+lo+m,real_temp);
    _mm256_store_pd(x->imag+lo+m,imag_temp);
    
  }
  if(n == 4)
  { 
    int m = 2;
    //We will first do the LEFT branch
    temp_real = x->real[lo];
    temp_im = x->imag[lo];
    x->real[lo] = temp_real + x->real[lo+1];
    x->imag[lo] = temp_im + x->imag[lo+1];
    x->real[lo+1] = temp_real - x->real[lo+1];
    x->imag[lo+1] = temp_im - x->imag[lo+1];
    //Now we do the RIGHT branch
    lo = lo+m;
    temp_real = x->real[lo];
    temp_im = x->imag[lo];
    //x->real[lo] = a-d, x->imag[lo] = b+c, x->real[lo+1] = a+d, x->imag[lo+1] = b-c
    //2a = (a-d)+(a+d), 2b = (b+c)+(b-c), 2c = (b+c)-(b-c), 2d = (a+d)-(a-d) 
    x->real[lo]   = temp_real     + x->real[lo+1];
    x->imag[lo]   = temp_im       + x->imag[lo+1];
    temp_real     = x->real[lo+1] - temp_real;
    x->real[lo+1] = temp_im       - x->imag[lo+1];
    x->imag[lo+1] = temp_real;
    //NOW WE WILL TIE EVERYTHING TOGETHER WE WILL GO FROM (x^2-1) and (x^2+1) TO (x^4-1)
    lo = lo-m;
    real_x_small = _mm_load_pd(x->real+lo);
    imag_x_small = _mm_load_pd(x->imag+lo);
    real_y_small = _mm_load_pd(x->real+lo+m);
    imag_y_small = _mm_load_pd(x->imag+lo+m);

    real_temp_small = _mm_sub_pd(real_x_small,real_y_small);
    imag_temp_small = _mm_sub_pd(imag_x_small,imag_y_small);
    real_x_small = _mm_add_pd(real_x_small,real_y_small);
    imag_x_small = _mm_add_pd(imag_x_small,imag_y_small);

    _mm_store_pd(x->real+lo,real_x_small);
    _mm_store_pd(x->imag+lo,imag_x_small);
    _mm_store_pd(x->real+lo+m,real_temp_small);
    _mm_store_pd(x->imag+lo+m,imag_temp_small);
  }
}