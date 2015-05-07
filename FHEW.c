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
    Ring_FFT ai;
    ct_ModQ res;
    int mm = (((m % q) + q) % q) * (2*N/q);             // Reduce mod q (dealing with negative number as well)
    int sign = 1;
    
    if (mm >= N) { 
    	mm -= N; 
    	sign = -1; 
    }
    
    for (int i = 0; i < K2; ++i) 
    {
      for (int k = 0; k < N; ++k) res[i][0][k] = rand(); // % Q
      	FFTforward(ai, res[i][0]);
      for (int k = 0; k < N2; ++k) 
        ai[k] = ((double complex) ai[k]) * ((double complex) sk_FFT[k]);
      FFTbackward(res[i][1], ai);
      for (int k = 0; k < N; ++k) 
      	res[i][1][k] += Sample(Chi1);    // Add error [a,as+e]
    }
    for (int i = 0; i < K; ++i) 
    {
      res[2*i  ][0][mm] += sign*vgprime[i]; // Add G Multiple
      res[2*i+1][1][mm] += sign*vgprime[i]; // [a,as+e] + X^m *G
    }
    for (int i = 0; i < K2; ++i)
      for (int j = 0; j < 2; ++j)
		FFTforward(ct[i][j], res[i][j]);

}


void KeyGen(EvalKey* EK, const SecretKey LWEsk) {
	SecretKeyN FHEWsk;
	KeyGenN(FHEWsk);
	SwitchingKeyGen(EK -> KSkey, LWEsk, FHEWsk);

	Ring_FFT FHEWskFFT;
	FFTforward(FHEWskFFT,FHEWsk);
	for (int i = 0; i < n; ++i)
      for (int j = 1; j < BS_base; ++j)
        for (int k = 0; k < BS_exp; ++k) 
        {
          EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));//IS THIS NEEDED?
          FHEWencrypt( (*EK->BSkey[i][j][k]), FHEWskFFT, LWEsk[i] * j * BS_table[k] );
        }

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
	        assert(fwrite(EK.KSkey[i][j][k], sizeof(CipherTextQ), 1, f));
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
		  EK->KSkey[i][j][k] = new CipherTextQ;
		  assert(fread(EK->KSkey[i][j][k], sizeof(CipherTextQ), 1, f));
		}
	return EK;
}

//TODOOOOOOOOOOOOOOOOOOO
// void AddToACC(ct_FFT ACC, ct_FFT C) {

// }

// void InitializeACC(ct_FFT ACC, int m) {

// }


// CipherTextQN* MemberTest(Ring_FFT t, ct_FFT C) {

// }

// void HomNAND(CipherText* res, const EvalKey& EK, const CipherText& ct1, const CipherText& ct2) {

// }

