#include <stdio.h>
#include "ring.h"
#include "fft.h"


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
