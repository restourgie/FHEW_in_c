#include <stdio.h>
#include "ring.h"
#include "fft.h"
#include <math.h>


static int ringdef_printed = 0;

static void ringdef_print()
{
  printf("R = IntegerModRing(2**32)\n");
  printf("RR = RR  = PolynomialRing(R,x)\n");
  printf("RRR = PolynomialQuotientRing(RR,RR(x**1024+1))\n");
  ringdef_printed = 1;
}

void ring_mul(ring_t *r, const ring_t *x, const ring_t *y)
{

  /* FFT based multiplication. */
  
  // double complex rfft[513];
  // double complex xfft[513];
  // double complex yfft[513];
  // int i;
  double complex rfft[2048];
  double complex xfft[2048];
  double complex yfft[2048];
  int i;

  FFTforward(xfft,x->v);
  FFTforward(yfft,y->v);

  for(i=0;i<513;i++)
    rfft[i] = xfft[i] * yfft[i];

  // for(i=0;i<2048;i++)
  //   rfft[i] = xfft[i] * yfft[i];

  FFTbackward(r->v,rfft);
}

void exact_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  /* Very simple schoolbook multiplication. Works. */
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

void to_complex(const ring_t *x, double complex *cplx_x)
{
  for(int i=0;i<CPLXDIM;++i)
    cplx_x[i] = x->v[i] + I*x->v[i+CPLXDIM];

  // for(int i=0;i<512;++i){
  //   printf("cplx_x[%d] = %f + I * %f\n x[%d] = %d x[%d] = %d\n",i,creal(cplx_x[i]),cimag(cplx_x[i]),i,x->v[i],i+512,x->v[i+512]);   
  // }
}

void to_real(const double complex *cplx_x, ring_t *x)
{
  for(int i=0;i<CPLXDIM;++i){
    x->v[i] = round(creal(cplx_x[i]));
    x->v[i+CPLXDIM] = round(cimag(cplx_x[i]));
  }
  // for(int i=0;i<512;++i){
  //     printf("cplx_x[%d] = %f + I * %f\n x[%d] = %u x[%d] = %u\n",i,creal(cplx_x[i]),cimag(cplx_x[i]),i,x->v[i],i+512,x->v[i+512]);
  // }
}

void naive_cplx_mul(double complex *res,const double complex *x, const double complex *y)
{
    double complex t;
    double complex big[REALDIM];

    for(int i=0;i<REALDIM;++i)
      big[i] = 0;

    for(int i =0;i<CPLXDIM;++i)
      for(int j=0;j<CPLXDIM;++j)
        big[i+j] += x[i] * y[j];    

    for(int i=CPLXDIM;i<REALDIM;++i){
      // printf("%f + i%f\n",creal(big[i+512]),cimag(big[i+512]));
      t = big[i] * (0+I*1);
      // printf("%f + i%f\n",creal(t),cimag(t));
      res[i-CPLXDIM] = big[i-CPLXDIM] + t;
    }
    
}

void fft_complex_mul(double complex *res, double complex *x,double complex *y)
{
  CPLX_FFTforward(x);
  CPLX_FFTforward(y);

  for(int i=0;i<CPLXDIM;++i)
    res[i] = x[i] * y[i];

  CPLX_FFTbackward(res);
  for(int i=0;i<CPLXDIM;++i)
    res[i] = res[i]/CPLXDIM;
}

void print_complex(const double complex *a){
    for(int i=0;i<CPLXDIM;i++)
      printf("cplxpoly[%d] = %f + i * %f\n",i,creal(a[i]),cimag(a[i]));
    printf("\n");
}

void complex_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_r[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);
  naive_cplx_mul(cplx_r,cplx_x,cplx_y);
  fft_complex_mul(cplx_res,cplx_x,cplx_y);
  printf("\n\nCOMPARING RESULTS\n\n");
  printf("Naive Poly\n");
  print_complex(cplx_r);
  printf("FFT complex\n");
  print_complex(cplx_res);
  to_real(cplx_r,r);
}


void ring_print(const ring_t *x)
{
  int i;
  if(!ringdef_printed)
    ringdef_print();

  printf("RRR(");
  for(i=1023;i>0;i--)
    printf("%u*x**%d + ",x->v[i],i);
  printf("%u*x**%d)",x->v[i],i);
}
