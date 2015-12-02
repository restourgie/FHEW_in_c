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
  init();
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t r,re,x,y;
  int n,i;
  int success = 0;

  for(n=0;n<NTESTS;n++)
  { 

    for(i=0;i<REALDIM;++i){
      x.v[i] = fgetc(urandom);
      y.v[i] = fgetc(urandom);
    }

     // x.v[0] = 58;
     // x.v[1] =37;
     // x.v[2] =238;
     // x.v[3] =155;
     // x.v[4] =1;
     // x.v[5] =153;
     // x.v[6] =239;
     // x.v[7] =69;
     // x.v[8] =239;
     // x.v[9] =152;
     // x.v[10] =60;
     // x.v[11] =43;
     // x.v[12] =148;
     // x.v[13] = 235;
     // x.v[14] = 149;
     // x.v[15] = 77;

     // y.v[0] = 58;
     // y.v[1] =37;
     // y.v[2] =238;
     // y.v[3] =155;
     // y.v[4] =1;
     // y.v[5] =153;
     // y.v[6] =239;
     // y.v[7] =69;
     // y.v[8] =239;
     // y.v[9] =152;
     // y.v[10] =60;
     // y.v[11] =43;
     // y.v[12] =148;
     // y.v[13] = 235;
     // y.v[14] = 149;
     // y.v[15] = 77;

    // fftw_mul(&re,&x,&y);
    // sr_vector_nonrec_mul(&re,&x,&y);
    // fftw_nega_mul(&re,&x,&y);
    negacyc_mul(&re,&x,&y);
    // sr_vector_mul(&re, &x, &y);
    naive_real_mul(&r,&x,&y);
    bool error = false;
    
    for(i=0;i<REALDIM;i++)
    {

      if(myabs(r.v[i] - re.v[i]) > THRESHOLD){
        // printf("school: %u \n", r.v[i]);
        // printf("mine: %u \n",re.v[i]);
        // printf("difference: %u\n\n",myabs(r.v[i] - re.v[i]));
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
        naive_real_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 1){
        start = rdtsc();
        negacyc_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 2){
        start = rdtsc();
        sr_vector_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 3){
        start = rdtsc();
        sr_vector_nonrec_mul(&r, &x, &y);
        end = rdtsc();
      }
      else if(j == 4){
        start = rdtsc();
        fftw_mul(&r,&x,&y);
        end = rdtsc();
      }
      else if(j == 5){
        start = rdtsc();
        fftw_nega_mul(&r,&x,&y);
        end = rdtsc();
      }
      cycles[n] = end - start;
    }
    qsort(cycles,sizeof(cycles)/sizeof(*cycles),sizeof(*cycles),compare);
    if(j == 0)
      printf("Naive real mul                     :%llu %llu %llu\n",cycles[249],cycles[499],cycles[749]);
    else if(j == 1)
      printf("Median Vectorized normal FFT       :%llu %llu %llu\n",cycles[249],cycles[499],cycles[749]);
    else if(j == 2)
      printf("Median Vectorized SplitRadix       :%llu %llu %llu\n",cycles[249],cycles[499],cycles[749]);
    else if(j == 3)
      printf("Median Vectorized Non recursive SR :%llu %llu %llu\n",cycles[249],cycles[499],cycles[749]);
    else if(j == 4)
      printf("Median FFTW                        :%llu %llu %llu\n",cycles[249],cycles[499],cycles[749]);
    else if(j == 5)
      printf("Median FFTW Dan's trick            :%llu %llu %llu\n",cycles[249],cycles[499],cycles[749]);
  }

  fclose(urandom); 
}

int main()
{
  // rand_test();
  cycle_meassure();

  return 0;
}