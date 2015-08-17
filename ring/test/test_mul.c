#include <stdio.h>
#include "../ring.h"

#define NTESTS 1000
#define THRESHOLD 2

void exact_mul(ring_t *r, const ring_t *x, const ring_t *y)
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
}

uint32_t myabs(uint32_t a)
{
  return (a>(1UL<<31)?-a:a);
}

int main()
{
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t r,re,x,y;
  int n,i;

  for(n=0;n<NTESTS;n++)
  { 
    fread(x.v,sizeof(uint32_t),1024,urandom);
    fread(y.v,sizeof(uint32_t),1024,urandom);
    for(i=0;i<1024;i++)
      y.v[i] &= 0x3ff;

    ring_mul(&r,&x,&y);
    exact_mul(&re,&x,&y);

    for(i=0;i<1024;i++)
    {
      if(myabs(r.v[i] - re.v[i]) > THRESHOLD)
        printf("%u\n",myabs(r.v[i] - re.v[i]));
    }
  }

  fclose(urandom);
  return 0;
}
