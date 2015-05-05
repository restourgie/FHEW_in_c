#include "FHEW.h"
#include "FFT.h"
#include <assert.h>

typedef Ring_ModQ ct_ModQ[K2][2];   // Ciphertext in coefficient form
typedef Ring_ModQ dct_ModQ[K2][K2]; // Decomposed Ciphertext in coeff form
typedef Ring_FFT  dct_FFT[K2][K2];  // Decomposed Ciphertext in FFT form

Ring_FFT t_TestMSB;

void setup(){
	Ring_modQ tmsb;
	for (int i = 1; i < N; ++i)
		tmsb[i]=1;
	FFTforward(t_TestMSB,tmsb);
}


void FHEWencrypt(ct_FFT ct, Ring_FFT sk_FFT, int m) {

}


void KeyGen(EvalKey* EK, const LWE::SecretKey LWEsk) {

}

void fwrite_ek(const EvalKey& EK, FILE* f) {
	// Write bootstrapping key
	for (int i = 0; i < n; ++i)      
	  for (int j = 1; j < BS_base; ++j)
	    for (int k = 0; k < BS_exp; ++k) 
	      assert(fwrite(EK.BSkey[i][j][k], sizeof(ct_FFT), 1, f));
	// Write switching key
	for (int i = 0; i < N; ++i)
	  for (int j = 0; j < KS_base; ++j)
	    for (int k = 0; k < KS_exp; ++k)
	        assert(fwrite(EK.KSkey[i][j][k], sizeof(LWE::CipherTextQ), 1, f));
}

EvalKey* fread_ek(FILE* f) {
	EvalKey* EK = new EvalKey;
	// Read bootstrapping key
	for (int i = 0; i < n; ++i)
	  for (int j = 1; j < BS_base; ++j)
	    for (int k = 0; k < BS_exp; ++k) 
	    {
	      //EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));
	      assert(fread(EK->BSkey[i][j][k], sizeof(ct_FFT), 1, f));
	    }
	// Read switching key
	printf("BSKey Read. \n");
	for (int i = 0; i < N; ++i)
	  for (int j = 0; j < KS_base; ++j)
	    for (int k = 0; k < KS_exp; ++k) 
	    {
		  EK->KSkey[i][j][k] = new LWE::CipherTextQ;
		  assert(fread(EK->KSkey[i][j][k], sizeof(LWE::CipherTextQ), 1, f));
		}
	return EK;
}


void AddToACC(ct_FFT ACC, ct_FFT C) {

}

void InitializeACC(ct_FFT ACC, int m) {

}


LWE::CipherTextQN* MemberTest(Ring_FFT t, ct_FFT C) {

}

void HomNAND(LWE::CipherText* res, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2) {

}

