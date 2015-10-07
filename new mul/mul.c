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

#define W(N,k) (cexp(2.0 * M_PI * I * (double)k / (double) N))
#define calc_cos(N,k) (cos(2.0 * M_PI * (double)k / (double) N))
#define calc_sin(N,k) (sin(2.0 * M_PI * (double)k / (double) N))

double **table;

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

// void print_cplx(cplx_ptr x,int n){
//   for (int i = 0; i < n; ++i)
//     printf("cplxpoly[%d] = %f + i * %f\n",i,x.real[i],x.imag[i]);
//   printf("\n");
// }

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

void twist(double complex *cplx_x,int n,int m,int lo)
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

void untwist(double complex *cplx_x,int n,int m,int lo)
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
* LOOKUPTABLE TWIST
*
******************************************************************/
void table_twist(cplx *cplx_x,int n,int m,int lo)
{
  double a,b,c,d;
  int j = 1, scale = ROOTDIM/n;
  for (int i = lo+1; i < lo+m; ++i)
  { 
    a = cplx_x->real[i];
    b = cplx_x->imag[i];  
    c = table[0][scale*j];
    d = table[1][scale*j];
     
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    cplx_x->real[i] = (a*c) - (b*d);
    cplx_x->imag[i] = (a*d) + (b*c);
    ++j;
  }
}

void table_untwist(cplx *cplx_x,int n,int m,int lo)
{
  double a,b,c,d;
  int j = 1, scale = ROOTDIM/n;
  for (int i = lo+1; i < lo+m; ++i)
  {
    a = cplx_x->real[i];
    b = cplx_x->imag[i];  
    c = table[0][scale*j];
    d = -table[1][scale*j];
     
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    cplx_x->real[i] = (a*c) - (b*d);
    cplx_x->imag[i] = (a*d) + (b*c);
    ++j;
  }
}

void init_table(){

  table = malloc(sizeof *table * 3);
  posix_memalign((void**)&table[0],64, ROOTDIM * sizeof(double));
  posix_memalign((void**)&table[1],64, ROOTDIM * sizeof(double));
  posix_memalign((void**)&table[2],64, ROOTDIM * sizeof(double));
  for (int i = 0; i < CPLXDIM; ++i)
  {
    table[0][i] = calc_cos(ROOTDIM,i);
    table[1][i] = calc_sin(ROOTDIM,i);
    table[2][i] = -table[1][i];
  }
}

/*
* This function initializes the lookup tables(LUTs) 
* We calculate all the roots of unity for 8,16,32,64,128,256 and 512
* We also calculate the first 512 roots of unity for 2048 
* This is needed for Bernsteins trick (x^512-i) equals (x^512-1) when twisted
* with the roots of 2048
*/
void init_table_vctr(){

  int size = 8, j=8;
  LUT1 = malloc(sizeof *LUT1 * size);
  LUT2 = malloc(sizeof *LUT2 * size);
  LUT3 = malloc(sizeof *LUT3 * size);
  for (int i = 0; i < size-1; ++i)
  {
    posix_memalign((void**)&LUT1[i],32, j * sizeof(double));
    posix_memalign((void**)&LUT2[i],32, j * sizeof(double));
    posix_memalign((void**)&LUT3[i],32, j * sizeof(double));
    for (int root = 0; root < j; ++root)
    {
      LUT1[i][root] = calc_cos(j,root);
      LUT2[i][root] = calc_sin(j,root);
      LUT3[i][root] = -LUT2[i][root];
    }
    j = j<<1;
  }
  posix_memalign((void**)&LUT1[7],32, 512 * sizeof(double));
  posix_memalign((void**)&LUT2[7],32, 512 * sizeof(double));
  posix_memalign((void**)&LUT3[7],32, 512 * sizeof(double));
  
  for (int i = 0; i < CPLXDIM; ++i)
  {
      LUT1[7][i] = calc_cos(ROOTDIM,i);
      LUT2[7][i] = calc_sin(ROOTDIM,i);
      LUT3[7][i] = -LUT2[7][i]; 
  }
}

/******************************************************************
*
* LOOKUPTABLE VECTOR TWIST
*
******************************************************************/
void vector_twist(cplx_ptr *cplx_x,int n,int m,int lo)
{
  __m256d real_x,imag_x,real_tbl,imag_tbl,imag_temp,real_temp;
  int j = 0, scale;
  if(n==ROOTDIM)
    scale = 7;
  else
    scale = log2(n)-3;
  for (int i = lo; i < lo+m; i+=4)
  { 
    // for (int bla = j; bla < j+4; ++bla)
    // {
    //   printf("root= %f + i * %f\n",table[0][scale*bla],table[1][scale*bla]);
    // }
    real_x = _mm256_load_pd(cplx_x->real+i);
    imag_x = _mm256_load_pd(cplx_x->imag+i);
    real_tbl = _mm256_load_pd(&LUT1[scale][j]);
    imag_tbl = _mm256_load_pd(&LUT2[scale][j]);
     
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_tbl);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_tbl);
    //REPLACED FOR COMMENTED SECTION
    //real_x = ac
    real_x = _mm256_mul_pd(real_x,real_tbl);
    //imag_x = bc
    imag_x = _mm256_mul_pd(imag_x,real_tbl);
    //real_x = ac - bd => real_x - real_temp
    real_x = _mm256_sub_pd(real_x,real_temp);
    //imag_x = ad + bc => imag_temp + imag_x
    imag_x = _mm256_add_pd(imag_x,imag_temp);
    //THESE ARE NOT WORKING 
    // real_x = _mm256_fmsub_pd(real_x,real_tbl,real_temp);
    // imag_x = _mm256_fmadd_pd(imag_x,real_tbl,imag_temp);
    _mm256_store_pd(cplx_x->real+i,real_x);
    _mm256_store_pd(cplx_x->imag+i,imag_x);
    j+=4;
  }
}

void vector_untwist(cplx_ptr *cplx_x,int n,int m,int lo)
{
  __m256d real_x,imag_x,real_tbl,imag_tbl,imag_temp,real_temp;
  int j = 0, scale;
  if(n==ROOTDIM)
    scale = 7;
  else
    scale = log2(n)-3;
  for (int i = lo; i < lo+m; i+=4)
  { 
    // for (int bla = j; bla < j+4; ++bla)
    // {
    //   printf("root= %f + i * %f\n",table[0][scale*bla],table[1][scale*bla]);
    // }
    real_x = _mm256_load_pd(cplx_x->real+i);
    imag_x = _mm256_load_pd(cplx_x->imag+i);
    real_tbl = _mm256_load_pd(&LUT1[scale][j]);
    imag_tbl = _mm256_load_pd(&LUT3[scale][j]);
     
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_tbl);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_tbl);
    //REPLACED FOR COMMENTED SECTION
    //real_x = ac
    real_x = _mm256_mul_pd(real_x,real_tbl);
    //imag_x = bc
    imag_x = _mm256_mul_pd(imag_x,real_tbl);

    //real_x = ac - bd => real_x - real_temp
    real_x = _mm256_sub_pd(real_x,real_temp);
    //imag_x = ad + bc => imag_temp + imag_x
    imag_x = _mm256_add_pd(imag_x,imag_temp);
    //THESE ARE NOT WORKING 
    // real_x = _mm256_fmsub_pd(real_x,real_tbl,real_temp);
    // imag_x = _mm256_fmadd_pd(imag_x,real_tbl,imag_temp);
    _mm256_store_pd(cplx_x->real+i,real_x);
    _mm256_store_pd(cplx_x->imag+i,imag_x);
    j+=4;
  }
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

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    cplx_x.real[i] = x->v[i];
    cplx_y.real[i] = y->v[i];
    cplx_x.imag[i] = x->v[j];
    cplx_y.imag[i] = y->v[j];
    ++j;
  }
  //TWIST CPLX_X AND APPLY FFT AFTERWARDS
  // vector_twist(&cplx_x,ROOTDIM,CPLXDIM,0);
  sr_vector(&cplx_x,CPLXDIM,0);
  // printf("\n\n**************VECTOR FFT X**************\n");
  // print_cplx(cplx_x,CPLXDIM);
  sr_vector_inverse(&cplx_x,CPLXDIM,0);
  // printf("\n\n**************VECTOR FFT X INVERSE**************\n");
  // print_cplx(cplx_x,CPLXDIM);
  //TWIST CPLX_Y AND APPLY FFT AFTERWARDS
  vector_twist(&cplx_y,ROOTDIM,CPLXDIM,0);
  sr_vector(&cplx_y,CPLXDIM,0);
  // print_cplx(cplx_y,CPLXDIM);

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
  sr_vector_inverse(&cplx_res,CPLXDIM,0);
  vector_untwist(&cplx_res,ROOTDIM,CPLXDIM,0);
  // printf("\n\n**************MULT RES**************\n");
  // print_cplx(cplx_res,CPLXDIM);

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
    // printf("\n\n**************X AFTER FFT**************\n");
  // print_double(&cplx_x,CPLXDIM);
  // sr_precomp(&cplx_x,CPLXDIM,0);
  // printf("\n\n**************X AFTER FFT**************\n");
  // print_double(&cplx_x,CPLXDIM);
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
  // printf("\n\n**************FFT NORMAL**************\n");	
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);
  // printf("\n\n**************Normal X**************\n");
  // print_complex(cplx_x,CPLXDIM);
  // twist(cplx_x,2*REALDIM,CPLXDIM,0);
  // printf("\n\n**************TWISTED X**************\n");
  // print_complex(cplx_x,CPLXDIM);
  twist(cplx_y,2*REALDIM,CPLXDIM,0);

  split_radix_recursive(cplx_x,CPLXDIM,0);
  // printf("\n\n**************FFT X**************\n");
  // print_complex(cplx_x,CPLXDIM);
  split_radix_recursive_inverse(cplx_x,CPLXDIM,0);
  // printf("\n\n**************FFT X INVERSE**************\n");
  // print_complex(cplx_x,CPLXDIM);

  split_radix_recursive(cplx_y,CPLXDIM,0);
  // printf("\n\n**************FFT X**************\n");
  // print_complex(cplx_y,CPLXDIM);

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
