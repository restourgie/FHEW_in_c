#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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

void print_double(const cplx *x,int N){
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
  table = malloc(sizeof *table * 2);
  table[0] = malloc(sizeof *table[0] *ROOTDIM);
  table[1] = malloc(sizeof *table[1] *ROOTDIM);
  for (int i = 0; i < CPLXDIM; ++i)
  {
    table[0][i] = calc_cos(ROOTDIM,i);
    table[1][i] = calc_sin(ROOTDIM,i);
  }
}

/******************************************************************
*
* LOOKUPTABLE VECTOR TWIST
*
******************************************************************/
void vector_twist(cplx *cplx_x,int n,int m,int lo)
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

void vector_untwist(cplx *cplx_x,int n,int m,int lo)
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
/******************************************************************
*
* SPLIT RADIX PRECOMPUTED AND VECTORIZED FFT MULTIPLICATION
*
******************************************************************/
void sr_vector_mul(ring_t *r, const ring_t *x, const ring_t *y){
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
  sr_vector(&cplx_x,CPLXDIM,0);
  
  // printf("\n\n**************X AFTER FFT**************\n");
  // print_double(&cplx_x,CPLXDIM);
  table_twist(&cplx_y,ROOTDIM,CPLXDIM,0);

  sr_vector(&cplx_y,CPLXDIM,0);

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
  table_untwist(&cplx_res,ROOTDIM,CPLXDIM,0);
  // // printf("\n\n**************MULT RES**************\n");
  // // print_double(&cplx_res,CPLXDIM);

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
  sr_precomp(&cplx_x,CPLXDIM,0);
  
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
  // // print_double(&cplx_res,CPLXDIM);

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
  twist(cplx_x,2*REALDIM,CPLXDIM,0);
  // printf("\n\n**************TWISTED X**************\n");
  // print_complex(cplx_x,CPLXDIM);
  twist(cplx_y,2*REALDIM,CPLXDIM,0);

  split_radix_recursive(cplx_x,CPLXDIM,0);
  // printf("\n\n**************FFT X**************\n");
  // print_complex(cplx_x,CPLXDIM);

  split_radix_recursive(cplx_y,CPLXDIM,0);

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
