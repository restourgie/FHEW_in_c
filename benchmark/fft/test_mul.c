#include <stdio.h>
#include <stdlib.h>
#include "../mul.h"
#include <stdbool.h>

#define NTESTS 1000
#define NALGO 2
#define THRESHOLD 2

// uint32_t myabs(uint32_t a)
// {
//   return (a>(1UL<<31)?-a:a);
// }

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

void cycle_meassure(){
  init();
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t r,x,y;
  int n,i;
  uint64_t start,end; 
  uint64_t cycles[NTESTS];

  for (int j = 0; j < NALGO; ++j)
  {
    for(n=0;n<NTESTS;n++)
    { 

      for(i=0;i<REALDIM;++i){
        x.v[i] = fgetc(urandom);
        y.v[i] = fgetc(urandom);
      }
      if(j == 0){
        start = rdtsc();
        fftw_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 1){
        start = rdtsc();
        fftw_nega_mul(&r,&x,&y);
        end = rdtsc();
      }
      cycles[n] = end - start;
    }
    qsort(cycles,sizeof(cycles)/sizeof(*cycles),sizeof(*cycles),compare);
    if(j == 0)
      printf("Naive real mul: %llu\n",cycles[NTESTS/2-1]);
    else if(j == 1)
      printf("Median Vectorized normal FFT: %llu\n",cycles[NTESTS/2-1]);

  }

  fclose(urandom); 
}

int main()
{
  // rand_test();
  cycle_meassure();

  return 0;
}