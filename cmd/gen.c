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
	// printf("sizeof BootstrappingKey: %lu\n",sizeof(BootstrappingKey));
	// printf("EK.BSkey: %p\n",EK.BSkey);
	if(EK.BSkey == NULL) {
		fprintf(stderr, "BAD BAD BAD!\n");
		return -1;
	}
	
	EK.KSkey = (SwitchingKey*) malloc(sizeof(SwitchingKey));
	
	if(EK.KSkey == NULL) {
		fprintf(stderr, "EVEN WORSE!\n");
		return -1;
	}
	// printf("Starting setup\n");
	Setup();
	//printf("Starting LWEKeyGen\n");
	LWEKeyGen(LWEsk); //APPROVED! THIS ONE WORKS
	

	// printf("Starting FHEWKeyGen\n");
	// printf("the size of EK = %lu\n",sizeof(EvalKey));
	
	// printf("%d ",((*EK.KSkey)[1024][0][0]).b);

	FHEWKeyGen(&EK,LWEsk);

	// printf("Printing the HUGE SwitchingKey\n");
	// for(int i =0; i < N;++i){
	// 	//printf("\n\n**********i is %d**********\n\n",i );
	// 	for(int j=0; j< KS_base;++j)
	// 		for(int k=0; k < KS_exp;++k){
	// 			// printf("%d\n",((*EK.KSkey)[i][j][k]).b);
	// 			for(int l=0;l < n;++l)	
	// 				printf("%d\n",((*EK.KSkey)[i][j][k]).a[l]);
	// 		}
	// }

	printf("Printing the huge BSkey\n");
	for(int i =0; i< n; ++i)
	  for(int j=1;j < BS_base;++j)
	    for(int k=0;k < BS_exp; ++k)
	      for(int l=0; l < K2; ++l)
	        for(int a =0; a <2; ++a)
	          for(int b =0; b < N2; ++b)
	           {
	                printf("%f\n",((*EK.BSkey)[i][j][k][l][a][b][0]));
	           }

	// printf("Saving EvalKey\n");
	// SaveEvalKey(&EK,ek_fn);
	// printf("Saving SecretKey\n");
	// SaveSecretKey(LWEsk,sk_fn);
	free(EK.KSkey);
	free(EK.BSkey);
}
