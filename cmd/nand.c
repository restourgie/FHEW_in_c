#include "../FHEW.h"
#include "common.h"
#include <stdlib.h>
#include <stdio.h>


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

  HomNAND(ct3, EK,ct1,ct2);

  SaveCipherText(ct3,oct_fn);
  free(EK->BSkey);
  free(EK->KSkey);
  free(EK);
  free(ct1);
  free(ct2);
  free(ct3);
}
