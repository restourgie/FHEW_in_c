#include "../FHEW.h"
#include "common.h"
#include <stdlib.h>
#include <stdio.h>

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

void help(char* cmd) {
  printf("\nusage: %s  EvalKeyFileName InCTFileName1 InCTFileName2 OutCTFileName  \n\n  Perform Homomorphic NAND computation.\n\n",cmd);
  exit(0);
}


int main(int argc, char *argv[]) {
  if (argc != 5) 
  	help(argv[0]);
  char* ek_fn = argv[1]; 
  char* ict1_fn = argv[2]; 
  char* ict2_fn = argv[3]; 
  char* oct_fn = argv[4]; 

  Setup();

  EvalKey* EK;

  EK = LoadEvalKey(ek_fn);

  CipherText *ct1,*ct2;

  ct1 = LoadCipherText(ict1_fn);
  ct2 = LoadCipherText(ict2_fn);

  CipherText *ct3 = malloc(sizeof(CipherText));
    uint64_t start,end; 
  uint64_t cycles[1000];    
  for (int count = 0; count < 1000; ++count)
  { 
    printf("count = %d\n",count );
    start = rdtsc(); 
    HomNAND(ct3, EK,ct1,ct2);
    end = rdtsc();
    cycles[count] = end - start;
  }
  qsort(cycles,sizeof(cycles)/sizeof(*cycles),sizeof(*cycles),compare);
  printf("Median: %llu\n",cycles[499]);

  SaveCipherText(ct3,oct_fn);
  free(EK->BSkey);
  free(EK->KSkey);
  free(EK);
  free(ct1);
  free(ct2);
  free(ct3);
}
