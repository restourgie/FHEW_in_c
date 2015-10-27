#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mul.h"
#include "negacyclic.h"
#include "fftw.h"

#define NTESTS 1000
#define NALGO 2
#define THRESHOLD 2

cplx_ptr vector_x;

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
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t x;
  int n,i;
  uint64_t start,end,start_2,end_2; 
  uint64_t cycles[NTESTS];
  uint64_t cycles_2[NTESTS];

  init_negacyc();
  posix_memalign((void**)&vector_x.real,32, CPLXDIM * sizeof(double));
  posix_memalign((void**)&vector_x.imag,32, CPLXDIM * sizeof(double));

  FFTsetup();
  double complex res[CPLXDIM];

  for (int j = 0; j < NALGO; ++j)
  {
    for(n=0;n<NTESTS;n++)
    { 

      for(i=0;i<REALDIM;++i)
        x.v[i] = fgetc(urandom);

      if(j == 0){
        start = rdtsc();
        phi_forward(&vector_x,&x);
        end = rdtsc();
        start_2 = rdtsc();
        phi_backward(&vector_x,&x);
        end_2 = rdtsc();
      }
      else if(j == 1){
        start = rdtsc();
        FFTWforward(res,&x);
        end = rdtsc();
        start_2 = rdtsc();
        FFTWbackward(&x,res);
        end_2 = rdtsc();
      }
      cycles[n] = end - start;
      cycles_2[n] = end_2 - start_2;
    }
    qsort(cycles,sizeof(cycles)/sizeof(*cycles),sizeof(*cycles),compare);
    qsort(cycles_2,sizeof(cycles_2)/sizeof(*cycles_2),sizeof(*cycles_2),compare);
    if(j == 0){
      printf("FFT Forward: %llu\n",cycles[NTESTS/2-1]);
      printf("FFT Backward: %llu\n",cycles_2[NTESTS/2-1]);
    }
    else if(j == 1){
      printf("FFTW Forward: %llu\n",cycles[NTESTS/2-1]);
      printf("FFTW Backward: %llu\n",cycles_2[NTESTS/2-1]);
    }

  }

  fclose(urandom); 
}

int main()
{
  cycle_meassure();

  return 0;
}