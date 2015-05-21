#include "../LWE.h"
#include "../FHEW.h"
#include "common.h"
#include <stdlib.h>
#include <stdio.h>



void help(char* cmd) {
  printf("\nusage: %s  Message SecretKeyFileName CipherTextFileName  \n\n   Encrypt the Message under some SecretKey and store it in a File.\n\n",cmd);
  exit(0);
}


int main(int argc, char *argv[]) {
  if (argc != 4) 
  	help(argv[0]);
  int message = atoi(argv[1]);
  char* sk_fn = argv[2]; 
  char* ct_fn = argv[3]; 

  if (!((message ==0)||(message ==1))){
      printf(" The message must be 0 or 1.\n");
  exit(0);
  }

  SecretKey *SK = LoadSecretKey(sk_fn);
  CipherText *ct = malloc(sizeof(CipherText));
  
  Encrypt(ct, *SK, message);  
  SaveCipherText(ct,ct_fn);
  
  free(SK);
  free(ct);
}


  