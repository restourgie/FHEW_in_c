#include <stdio.h>
#include <stdlib.h>
#include "../mul.h"
#include <stdbool.h>

#define NTESTS 1000
#define THRESHOLD 2

uint64_t start,end; 
uint64_t cycles[NTESTS];

uint32_t myabs(uint32_t a)
{
  return (a>(1UL<<31)?-a:a);
}

int compare(const void * elem1, const void * elem2)
{	
	uint64_t x = *((uint64_t*)elem1);
	uint64_t y = *((uint64_t*)elem2);
	if(x > y) 
		return 1;
	else if(x < y) 
		return -1;
	return 0;
}

uint64_t rdtsc()
{
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
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
    // start = rdtsc();
    // normal_fft_mul(&r,&x,&y);
    // split_radix_mul(&re,&x,&y);
    sr_vector_mul(&r,&x,&y);
    fftw_mul(&re,&x,&y);
    // sr_precomp_mul(&re,&x,&y);
    // naive_complex_mul(&r,&x,&y);
    // naive_real_mul(&r,&x,&y);
    // end = rdtsc();
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
    // cycles[n] = end - start;
  }
  // printf("Ammount of successful multiplications: %d\n", success);

  fclose(urandom);
  // qsort(cycles,sizeof(cycles)/sizeof(*cycles),sizeof(*cycles),compare);
  // for (int i = 0; i < NTESTS; ++i)
  // {
	 //  printf("%lu\n",cycles[i]); //<< "index[" << i << "] : "
  // }
  // printf("middle: %lu\n",cycles[499]);

}

int main()
{
  rand_test();

  return 0;
}