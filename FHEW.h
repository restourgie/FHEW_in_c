#ifndef FHEW_H
#define FHEW_H

#include "params.h"
#include "LWE.h"
#include <stdio.h>

  typedef Ring_FFT ct_FFT[K2][2];  // Ciphertext in FFT form ct_FFT[K2][2] => Ring_FFT[6][2] => complex double[6][2][513]
  typedef ct_FFT BootstrappingKey[n][BS_base][BS_exp]; 
  //BootstrappingKey[500][23][2] => ct_FFT[500][23][2] => ct_FFT[500][23][2][6][2] => complex double[500][23][2][6][2][513]
  typedef struct {
    BootstrappingKey *BSkey;
    SwitchingKey *KSkey;
  } EvalKey;

  void Setup();
  void FHEWKeyGen(EvalKey* EK,SecretKey LWEsk);
  void HomNAND(CipherText* res, EvalKey* EK, CipherText* ct1, CipherText* ct2); 

#endif
