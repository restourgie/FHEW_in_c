#ifndef FHEW_H
#define FHEW_H

#include "params.h"
#include "LWE.h"
#include <stdio.h>

/*
	Ring_FFT = 513*2*double;
	ct_FFT = 6*2*Ring_FFT;
	BootstrappingKey = 500*23*2*ct_FFT;
*/

	// typedef struct {
	// 	Ring_FFT v[K2][2];
	// } ct_FFT;

  typedef Ring_FFT ct_FFT[K2][2];  // Ciphertext in FFT form ct_FFT[K2][2] => Ring_FFT[6][2] => complex_double[6][2][513] => double[6][2][513][2]
  typedef ct_FFT BootstrappingKey[n][BS_base][BS_exp]; 
  //BootstrappingKey[500][23][2] => ct_FFT[500][23][2] => ct_FFT[500][23][2][6][2] => complex_double[500][23][2][6][2][513] => double[500][23][2][6][2][513][2]
  typedef struct {
    BootstrappingKey *BSkey;
    SwitchingKey *KSkey;
  } EvalKey;

  void Setup();
  void FHEWKeyGen(EvalKey* EK,SecretKey LWEsk);
  void HomNAND(CipherText* res, EvalKey* EK, CipherText* ct1, CipherText* ct2); 

#endif
