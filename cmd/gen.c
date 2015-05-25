#include <stdio.h>
#include "../LWE.h"
#include "../FHEW.h"
#include "common.h"
#include <stdlib.h>

void help(char* cmd){
	printf("\nusage: %s  SecretKeyFileName EvalKeyFileName  \n\n Generate a secret key sk and evaluation key ek, and store them in two separate files.\n\n",cmd);
	exit(0);
}


int main (int argc, char *argv[])
{
	if(argc != 3) 
		help(argv[0]);

	char* sk_fn = argv[1]; 
  char* ek_fn = argv[2];

  SecretKey LWEsk;
  EvalKey *EK;
  EK = (EvalKey*) malloc(sizeof(EvalKey));
  printf("\n Pointer value of EK = %p\n",EK);
  printf("\n Pointer value of EK->BSkey = %p\n",&(EK->BSkey));
  printf("\n Pointer value of EK->KSkey = %p\n",&(EK->KSkey));

  printf("Starting setup\n");
  Setup();
  printf("Starting LWEKeyGen\n");
  LWEKeyGen(LWEsk);
  printf("Starting FHEWKeyGen\n");
  printf("the size of EK = %lu\n",sizeof(EvalKey));
  FHEWKeyGen(EK,LWEsk);
  printf("\n Pointer value of EK = %p\n",EK);
  printf("Saving EvalKey\n");
  // SaveEvalKey(EK,ek_fn);
  // printf("Saving SecretKey\n");
  // SaveSecretKey(LWEsk,sk_fn);
  free(EK);
}