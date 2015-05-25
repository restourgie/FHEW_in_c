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

void SaveSecretKey(SecretKey LWEsk, char* filepath) {
  FILE * f;
  f = fopen(filepath, "wb"); // wb -write binary
  if (f == NULL) {
    printf("Failed to open %s in Write-Binary mode .\n",filepath);
    exit(EXIT_FAILURE);
  }
  printf("Writing Secret key to %s .\n", filepath);
  fwrite(&LWEsk, sizeof(SecretKey), 1, f);
  fclose(f);
}

SecretKey* LoadSecretKey(char* filepath) {
  FILE * f;
  f = fopen(filepath, "rb"); // wb -write binary
  if (f == NULL) {
    printf("Failed to open %s in Read-Binary mode.\n", filepath);
    exit(EXIT_FAILURE);
  }
  SecretKey *LWEsk = malloc(sizeof(SecretKey));  
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
  
//    Write bootstrapping key
  for (int i = 0; i < n; ++i)      
    for (int j = 1; j < BS_base; ++j)
      for (int k = 0; k < BS_exp; ++k){ 

           //fwrite(&(EK->BSkey[i][j][k]), sizeof(ct_FFT), 1, f);
      }
    // for (int i = 0; i < n; ++i){
    // printf("i: %d\n",i);
    // for (int j = 1; j < BS_base; ++j){
    //   printf("j: %d\n",j );
    //   for (int k = 0; k < BS_exp; ++k){
    //     printf("k: %d\n", k);
    //     for(int l = 0; l < K2; ++l){
    //       printf("l: %d\n", l);
    //       for(int bin = 0; bin < 2; ++bin){
    //         printf("bin: %d\n",bin);
    //         for(int last = 0; last < N2; ++last){
    //           printf("last: %d\n",last);
    //           printf("BSkey real: %f imag: %f\n",EK->BSkey[i][j][k][l][bin][last][0],EK->BSkey[i][j][k][l][bin][last][1]);
    //         }}}}}}

  // Write switching key
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < KS_base; ++j)
      for (int k = 0; k < KS_exp; ++k)
          assert(fwrite(&(EK->KSkey[i][j][k]), sizeof(CipherTextQ), 1, f));

  fclose(f);
}

EvalKey* LoadEvalKey(char* filepath) {
  EvalKey *EK = malloc(sizeof(EvalKey));
  FILE * f;
  f = fopen(filepath, "rb"); // rb -read binary
  if (f == NULL){
    printf("Failed to open %s in Read-Binary mode .\n",filepath);
    exit(EXIT_FAILURE);
  }
  printf("Reading Evaluation key from %s.\n", filepath);

  // Read bootstrapping key
  for (int i = 0; i < n; ++i)
    for (int j = 1; j < BS_base; ++j)
      for (int k = 0; k < BS_exp; ++k) 
      {
        //EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));
        assert(fread(&(EK->BSkey[i][j][k]), sizeof(ct_FFT), 1, f));
      }
  printf("BSKey Read. \n");
  
  // Read switching key
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < KS_base; ++j)
      for (int k = 0; k < KS_exp; ++k) 
      {
        // EK->KSkey[i][j][k] = new CipherTextQ;//??????????????????????????????????????????????????????????????
        assert(fread(&(EK->KSkey[i][j][k]), sizeof(CipherTextQ), 1, f));
      }
  printf("KSKey Read : %d \t %d \t %d .\n", N, KS_base, KS_exp);
  
  fclose(f);
  return EK;
}

/*************************************************************************
*                                                                        *
*                           CIPHER TEXT                                  *
*                                                                        *
*************************************************************************/


void SaveCipherText(CipherText* ct, char* filepath){
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
  CipherText *ct = malloc(sizeof(CipherText));
  assert(fread(ct, sizeof(CipherText), 1, f));
  fclose(f);
  return ct;
}