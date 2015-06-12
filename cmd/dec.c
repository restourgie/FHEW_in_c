#include "../FHEW.h"
#include "../LWE.h"
#include "common.h"
#include <stdlib.h>
#include <stdio.h>


void help(char* cmd) {
  printf("\nusage: %s  SecretKeyFileName CipherTextFileName  \n\n   Decrypt the CipherText under some SecretKey and print it on the std output.\n\n",cmd);
  exit(0);
}


int main(int argc, char *argv[]) {
	  if (argc != 3) 
      help(argv[0]);
  char* sk_fn = argv[1]; 
  char* ct_fn = argv[2]; 

  // Setup(); IS THIS NEEDED?????

  SecretKey* SK = LoadSecretKey(sk_fn);
  CipherText* ct = LoadCipherText(ct_fn);
  int m = Decrypt(*SK,ct);
  printf("%d\n",m);
  free(SK);
  free(ct);
}
