#ifndef FHEW_H
#define FHEW_H

#include "params.h"
#include "LWE.h"


  typedef Ring_FFT  ct_FFT[K2][2];  // Ciphertext in FFT form ct_FFT[K2][2] => Ring_FFT[6][2] => complex_double[6][2][513] => double[6][2][513][2]
  typedef ct_FFT* BootstrappingKey[n][BS_base][BS_exp]; 
  //BootstrappingKey[500][23][2] => ct_FFT[500][23][2] => ct_FFT[500][23][2][6][2] => complex_double[500][23][2][6][2][513] => double[500][23][2][6][2][513][2]
  typedef struct {
    BootstrappingKey BSkey;
    SwitchingKey KSkey;
  } EvalKey;

  void Setup();
  void FHEWKeyGen(EvalKey* EK, const SecretKey LWEsk);
  //TODOOOOOOOOOOOOOOOOO
  // void HomNAND(CipherText* res, const EvalKey& EK, const CipherText& ct1, const CipherText& ct2);

  //void fwrite_ek(const EvalKey& EK, FILE* f);
  // EvalKey* fread_ek(FILE* f);

#endif
