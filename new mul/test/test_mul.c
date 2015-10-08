#include <stdio.h>
#include "../mul.h"
#include <stdbool.h>

#define NTESTS 1000
#define THRESHOLD 2

uint32_t myabs(uint32_t a)
{
  return (a>(1UL<<31)?-a:a);
}

void rand_test(){
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t r,re,x,y;
  int n,i;
  int success = 0;
  init();

  for(n=0;n<NTESTS;n++)
  { 

    for(i=0;i<REALDIM;++i){
      x.v[i] = fgetc(urandom);
      y.v[i] = fgetc(urandom);
    }

    // normal_fft_mul(&r,&x,&y);
    // split_radix_mul(&re,&x,&y);
    sr_vector_mul(&r,&x,&y);
    fftw_mul(&re,&x,&y);
    // sr_precomp_mul(&re,&x,&y);
    // naive_complex_mul(&r,&x,&y);
    // naive_real_mul(&r,&x,&y);
    bool error = false;
    
    for(i=0;i<REALDIM;i++)
    {

      if(myabs(r.v[i] - re.v[i]) > THRESHOLD){
        printf("school: %u \n", re.v[i]);
        printf("mine: %u \n",r.v[i]);
        printf("difference: %u\n\n",myabs(r.v[i] - re.v[i]));
        error = true;
      }
      
    }
    if(!error)
      success++;

  }
  printf("Ammount of successful multiplications: %d\n", success);

  fclose(urandom);

}

int main()
{
  rand_test();

  return 0;
}