#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "mul.h"
#include "fft/fft_negacyc.h"
#include "fft/split_radix_fft.h"
#include "fft/sr_precomp.h"
#include "fft/sr_vector.h"
#include "fft/fftw.h"



/******************************************************************
*
* SUPPORT CODE
*
******************************************************************/
void print_complex(const double complex *a, int N){
    for(int i=0;i<N;++i)
      printf("cplxpoly[%d] = %f + i * %f\n",i,creal(a[i]),cimag(a[i]));
    printf("\n");
}

void print_double(const cplx *x,int N){
  for (int i = 0; i < N; ++i)
    printf("cplxpoly[%d] = %f + i * %f\n",i,x->real[i],x->imag[i]);
  printf("\n");
}

void print_cplx(const cplx_ptr *x,int N){
  for (int i = 0; i < N; ++i)
    printf("cplxpoly[%d] = %f + i * %f\n",i,x->real[i],x->imag[i]);
  printf("\n");
}

/******************************************************************
*
* CONVERSION
*
******************************************************************/
void to_complex(const ring_t *x, double complex *cplx_x)
{
  for(int i=0;i<CPLXDIM;++i)
    cplx_x[i] = x->v[i] + I*x->v[i+CPLXDIM];
}

void to_real(const double complex *cplx_x, ring_t *x)
{
  for(int i=0;i<CPLXDIM;++i){
    x->v[i] = round(creal(cplx_x[i]));
    x->v[i+CPLXDIM] = round(cimag(cplx_x[i]));
  }
}

void init(){
  //PRECOMP TABLES FOR VECTOR FFT
  init_table_vctr();
  //PRECOMP TABLES FOR PRECOMP FFT
  init_table();
  //PRECOMP FFTW
  FFTsetup();
}

/******************************************************************
*
* FFTW MULTIPLICATION
*
******************************************************************/
void fftw_mul(ring_t *r, const ring_t *x, const ring_t *y){
  double complex cplx_x[CPLXDIM+1];
  double complex cplx_y[CPLXDIM+1];
  double complex cplx_res[CPLXDIM+1];

  FFTWforward(cplx_x,x);
  FFTWforward(cplx_y,y);

  for (int i = 0; i < CPLXDIM+1; ++i)
  {
    cplx_res[i] = cplx_x[i] * cplx_y[i];
  }

  FFTWbackward(r,cplx_res);

}

/******************************************************************
*
* SPLIT RADIX PRECOMPUTED AND VECTORIZED FFT MULTIPLICATION
*
******************************************************************/
void sr_vector_mul(ring_t *r, const ring_t *x, const ring_t *y){
  // printf("\n\n**************split-radix FAST**************\n");
  cplx_ptr cplx_x,cplx_y,cplx_res;
  posix_memalign((void**)&cplx_x.real,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&cplx_x.imag,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&cplx_y.real,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&cplx_y.imag,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&cplx_res.real,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&cplx_res.imag,32, CPLXDIM * sizeof(double));


  fft_vector_forward(&cplx_x,x);
  // print_cplx(&cplx_x,CPLXDIM);
  fft_vector_forward(&cplx_y,y);
  
  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp;
  // double a,b,c,d;
  for (int i = 0; i < CPLXDIM; i+=4)
  {
    real_x = _mm256_load_pd(cplx_x.real+i);
    imag_x = _mm256_load_pd(cplx_x.imag+i);
    real_y = _mm256_load_pd(cplx_y.real+i);
    imag_y = _mm256_load_pd(cplx_y.imag+i);

    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_y);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_y);
    //REPLACED FOR COMMENTED SECTION
    //real_x = ac
    real_x = _mm256_mul_pd(real_x,real_y);
    //imag_x = bc
    imag_x = _mm256_mul_pd(imag_x,real_y);
    //real_x = ac - bd => real_x - real_temp
    real_x = _mm256_sub_pd(real_x,real_temp);
    //imag_x = ad + bc => imag_temp + imag_x
    imag_x = _mm256_add_pd(imag_x,imag_temp);
    //THESE ARE NOT WORKING 
    // real_x = _mm256_fmsub_pd(real_x,real_tbl,real_temp);
    // imag_x = _mm256_fmadd_pd(imag_x,real_tbl,imag_temp);
    _mm256_store_pd(cplx_res.real+i,real_x);
    _mm256_store_pd(cplx_res.imag+i,imag_x);

    // a = cplx_x.real[i];
    // b = cplx_x.imag[i]; 
    // c = cplx_y.real[i];
    // d = cplx_y.imag[i];
    // //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    // cplx_res.real[i] = ((a*c) - (b*d))/CPLXDIM;
    // cplx_res.imag[i] = ((a*d) + (b*c))/CPLXDIM;
  }
  // print_cplx(&cplx_res,CPLXDIM);
  fft_vector_backward(&cplx_res,r);
  
}

/******************************************************************
*
* SPLIT RADIX PRECOMPUTED FFT MULTIPLICATION
*
******************************************************************/
void sr_precomp_mul(ring_t *r, const ring_t *x, const ring_t *y){
  // printf("\n\n**************split-radix FAST**************\n");
  cplx cplx_x,cplx_y,cplx_res;

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x.real[i] = x->v[i];
    cplx_y.real[i] = y->v[i];
    cplx_x.imag[i] = x->v[j];
    cplx_y.imag[i] = y->v[j];
    ++j;
  }
  table_twist(&cplx_x,ROOTDIM,CPLXDIM,0);
  sr_precomp(&cplx_x,CPLXDIM,0);

  table_twist(&cplx_y,ROOTDIM,CPLXDIM,0);
  sr_precomp(&cplx_y,CPLXDIM,0);

  double a,b,c,d;
  for (int i = 0; i < CPLXDIM; ++i)
  {
      a = cplx_x.real[i];
      b = cplx_x.imag[i]; 
      c = cplx_y.real[i];
      d = cplx_y.imag[i];
      //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
      cplx_res.real[i] = ((a*c) - (b*d))/CPLXDIM;
      cplx_res.imag[i] = ((a*d) + (b*c))/CPLXDIM;
  }
  // // printf("\n\n**************STARTING INVERSE**************\n");
  sr_precomp_inverse(&cplx_res,CPLXDIM,0);
  table_untwist(&cplx_res,ROOTDIM,CPLXDIM,0);
  // // printf("\n\n**************MULT RES**************\n");
  // print_double(&cplx_res,CPLXDIM);

  j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    r->v[i] = cplx_res.real[i];
    r->v[j] = cplx_res.imag[i];
    ++j; 
  }
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

  
  twist(cplx_x,ROOTDIM,CPLXDIM,0);
  split_radix_recursive(cplx_x,CPLXDIM,0);

  twist(cplx_y,ROOTDIM,CPLXDIM,0);
  split_radix_recursive(cplx_y,CPLXDIM,0);


  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }
  split_radix_recursive_inverse(cplx_res,CPLXDIM,0);
  untwist(cplx_res,ROOTDIM,CPLXDIM,0);

  to_real(cplx_res,r);
}

/******************************************************************
*
* SCHOOLBOOK NEGACYCLIC FFT
*
******************************************************************/
void normal_fft_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);

  double complex root = I;
  root = csqrt(root);
  
  recursive_phi(cplx_x,CPLXDIM,0,root);
  recursive_phi(cplx_y,CPLXDIM,0,root);


  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  }

  inverse_phi(cplx_res,CPLXDIM,0,root);
  // printf("\n\n**************SMART COMPLEX MUL RESULT**************\n");
  // print_complex(cplx_res,CPLXDIM);
  to_real(cplx_res,r);
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
