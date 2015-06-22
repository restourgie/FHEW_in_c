#include <stdio.h>
#include "ring.h"


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
