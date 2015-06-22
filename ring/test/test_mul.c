#include <stdio.h>
#include "../ring.h"

#define NTESTS 10

int main()
{
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t r,x,y;
  int n;

  for(n=0;n<NTESTS;n++)
  {
    fread(x.v,sizeof(uint32_t),1024,urandom);
    fread(y.v,sizeof(uint32_t),1024,urandom);

    ring_mul(&r,&x,&y);

    ring_print(&x);
    printf(" * ");
    ring_print(&y);
    printf(" - ");
    ring_print(&r);
    printf("\n");
  }

  fclose(urandom);
  return 0;
}
