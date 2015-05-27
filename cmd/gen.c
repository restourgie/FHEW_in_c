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
	EvalKey EK;
	EK.BSkey = (BootstrappingKey*) malloc(sizeof(BootstrappingKey));
	printf("sizeof BootstrappingKey: %lu\n",sizeof(BootstrappingKey));
	printf("EK.BSkey: %p\n",EK.BSkey);
	if(EK.BSkey == NULL) {
		fprintf(stderr, "BAD BAD BAD!\n");
		return -1;
	}
	/* TODO: check return value of malloc */
	EK.KSkey = (SwitchingKey*) malloc(sizeof(SwitchingKey));
	/* TODO: check return value of malloc */
	if(EK.KSkey == NULL) {
		fprintf(stderr, "EVEN WORSE!\n");
		return -1;
	}
	printf("\nPointer value of EK.BSkey = %p\n",(EK.BSkey));
	printf("\nPointer value of EK.KSkey = %p\n",(EK.KSkey));

	printf("Starting setup\n");
	Setup();
	printf("Starting LWEKeyGen\n");
	LWEKeyGen(LWEsk);
	printf("Starting FHEWKeyGen\n");
	printf("the size of EK = %lu\n",sizeof(EvalKey));
	FHEWKeyGen(&EK,LWEsk);
	printf("Saving EvalKey\n");
	// SaveEvalKey(EK,ek_fn);
	// printf("Saving SecretKey\n");
	// SaveSecretKey(LWEsk,sk_fn);
	free(EK.KSkey);
	free(EK.BSkey);
}
