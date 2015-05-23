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

  EvalKey *EK;
  SecretKey LWEsk;

  printf("Starting setup\n");
  Setup();
  printf("Starting LWEKeyGen\n");
  LWEKeyGen(LWEsk);
  printf("Starting FHEWKeyGen\n");
  EK = FHEWKeyGen(LWEsk);
  printf("Saving EvalKey\n");
  SaveEvalKey(EK,ek_fn);
  printf("Saving SecretKey\n");
  SaveSecretKey(LWEsk,sk_fn);
  free(EK);
}