#include "../LWE.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../FHEW.h"


/*************************************************************************
*                                                                        *
*                           SECRET KEY                                   *
*                                                                        *
*************************************************************************/

void SaveSecretKey(const SecretKey* LWEsk, char* filepath) {
  FILE * f;
  f = fopen(filepath, "wb"); // wb -write binary
  if (f == NULL) {
    printf("Failed to open %s in Write-Binary mode .\n",filepath);
    exit(EXIT_FAILURE);
  }
  printf("Writing Secret key to %s .\n", filepath);
  fwrite(LWEsk, sizeof(SecretKey), 1, f);
  fclose(f);
}

SecretKey* LoadSecretKey(char* filepath) {
  FILE * f;
  f = fopen(filepath, "rb"); // wb -write binary
  if (f == NULL) {
    printf("Failed to open %s in Read-Binary mode.\n", filepath);
    exit(EXIT_FAILURE);
  }
  SecretKey* LWEsk = (SecretKey*) malloc(sizeof(SecretKey));  
  printf("Reading Secret key From %s .\n", filepath);
  assert(fread(LWEsk, sizeof(SecretKey), 1, f));
  printf("Secret Key read.\n");
  fclose(f);
  return LWEsk;
 
}

/*************************************************************************
*                                                                        *
*                           EVAL KEY                                     *
*                                                                        *
*************************************************************************/

void SaveEvalKey(EvalKey *EK, char* filepath) {
  FILE * f;
  f = fopen(filepath, "wb"); // wb -write binary
  if (f == NULL) {
    printf("Failed to open %s in Write-Binary mode .\n", filepath);
    exit(EXIT_FAILURE);
  }
  printf("Writing Evaluation key to %s .\n", filepath);
  fwrite_ek(*EK, f);
  fclose(f);
}

EvalKey* LoadEvalKey(char* filepath) {
  EvalKey* EK;
  FILE * f;
  f = fopen(filepath, "rb"); // rb -read binary
  if (f == NULL){
    printf("Failed to open %s in Read-Binary mode .\n",filepath);
    exit(EXIT_FAILURE);
  }
  printf("Reading Evaluation key from %s.\n", filepath);
  EK = fread_ek(f);
  printf("KSKey Read : %d \t %d \t %d .\n", N, KS_base, KS_exp);
  fclose(f);
  return EK;
}

/*************************************************************************
*                                                                        *
*                           CIPHER TEXT                                  *
*                                                                        *
*************************************************************************/


void SaveCipherText(const CipherText* ct, char* filepath){
   FILE * f;
  f = fopen(filepath, "wb"); // wb -write binary
  if (f == NULL){
    printf("Failed to open %s in Write-Binary mode .\n", filepath);
    exit(EXIT_FAILURE);
  }
  printf("Writing CipherText to %s .\n", filepath);
  assert(fwrite(ct, sizeof(CipherText), 1, f));
  fclose(f);
}

CipherText* LoadCipherText(char* filepath) {
  FILE * f;
  f = fopen(filepath, "rb"); // wb -write binary
  if (f == NULL) {
    printf("Failed to open %s in Write-Binary mode .\n", filepath);
    exit(EXIT_FAILURE);
  }
  printf("Loading CipherText from %s.\n", filepath);
  CipherText *ct = malloc(sizeof *ct);
  assert(fread(ct, sizeof(CipherText), 1, f));
  fclose(f);
  return ct;
}