#include <stdio.h>
#include <stdlib.h>
#include "../mul.h"
#include <stdbool.h>

#define NTESTS 1000
#define NALGO 6
#define THRESHOLD 2

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
    // normal_fft_mul(&re,&x,&y);
    // split_radix_mul(&re,&x,&y);
    sr_vector_mul(&r,&x,&y);
    sr_vector_nonrec_mul(&re,&x,&y);
    // fftw_mul(&re,&x,&y);
    // sr_precomp_mul(&re,&x,&y);
    // naive_complex_mul(&r,&x,&y);
    // naive_real_mul(&r,&x,&y);
    // test(&re,&x,&y);
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

void cycle_meassure(){
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t r,x,y;
  int n,i;
  uint64_t start,end; 
  uint64_t cycles[NTESTS];
  init();

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
        normal_fft_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 1){
        start = rdtsc();
        split_radix_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 2){
        start = rdtsc();
        sr_precomp_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 3){
        start = rdtsc();
        sr_vector_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 4){
        start = rdtsc();
        fftw_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 5){
        start = rdtsc();
        sr_vector_nonrec_mul(&r,&x,&y);
        end = rdtsc();
      }
      cycles[n] = end - start;
    }
    qsort(cycles,sizeof(cycles)/sizeof(*cycles),sizeof(*cycles),compare);
    if(j == 0)
      printf("Median Normal FFT: %llu\n",cycles[499]);
    else if(j == 1)
      printf("Median Split Radix FFT: %llu\n",cycles[499]);
    else if(j == 2)
      printf("Median Split Radix Precomputed: %llu\n",cycles[499]);
    else if(j == 3)
      printf("Median Split Radix VECTORIZED: %llu\n",cycles[499]);
    else if(j == 4)
      printf("Median FFTW: %llu\n",cycles[499]);
    else if(j == 5)
      printf("Median Split Radix VECTORIZED NON RECURSIVE: %llu\n",cycles[499]);
  }

  fclose(urandom); 
}

int main()
{
  // rand_test();
  cycle_meassure();

  return 0;
}