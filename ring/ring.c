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
  /* Very simple schoolbook multiplication. Works. */
  int i,j;
  for(i=0;i<1024;i++)
    r->v[i] = 0;
  for(i=0;i<1024;i++)
  {
    for(j=0;j<1024;j++)
    {
      if(i+j < 1024)
        r->v[i+j] += x->v[i] * y->v[j];
      else
        r->v[i+j-1024] -= x->v[i] * y->v[j];
    }
  }
  /* FFT based multiplication. Doesn't work. */
  /*
  double complex rfft[513];
  double complex xfft[513];
  double complex yfft[513];
  int i;

  FFTforward(xfft,x->v);
  FFTforward(yfft,y->v);

  for(i=0;i<513;i++)
    rfft[i] = xfft[i] * yfft[i];

  FFTbackward(r->v,rfft);
  */
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
