#include "FHEW.h"
#include "params.h"
#include "stdbool.h"
#include <stdlib.h>


void LWEKeyGen(SecretKey sk) {
	KeyGenRestart:;
    int s=0, ss=0;
    for (int i = 0; i < n; ++i) {

          sk[i] = Sample(Chi_Binary);
          s+= sk[i];
          ss+= abs(sk[i]);
        }
      if (abs(s)>5) goto KeyGenRestart;
      if (abs(ss - n/2)>5) goto KeyGenRestart;
}

void KeyGenN(SecretKeyN sk) {
	for(int i =0;i < N; ++i)
		sk[i] = Sample(Chi1); //WHY ARE THERE NO CHECKS HERE?? IT IS THE SAME
}

void Encrypt(CipherText* ct, const SecretKey sk, int m) {
    ct->b = (m % 4) * q / 4;//can you just do m * q / 4 because m is 0 or 1 this is always the same mod 4
    ct->b += Sample(Chi3);
    
    for (int i = 0; i < n; ++i)	
    {
      ct->a[i] = rand() % q;
      ct->b = (ct->b + ct->a[i] * sk[i]) % q;
    }
}

// int Decrypt(const SecretKey sk, const CipherText& ct) {
 
// }

// void DecryptDetail(const SecretKey sk, const CipherText& ct) {

// }

    
// int round_qQ(ZmodQ v) {
  
// }
  
// void ModSwitch(CipherText* ct, const CipherTextQ& c) {

// }

//GENERATE SwitchingKey //SwitchingKey => CipherTextQ[1024][25][7] //CipherTextQ = {ZmodQ a[n]; ZmodQ b;} => ZmodQ = int32_t and n = 500 
void SwitchingKeyGen(SwitchingKey res,const SecretKey new_sk,const SecretKeyN old_sk) {
	for (int i = 0; i < N; ++i) 
      for (int j = 0; j < KS_base; ++j)
		for (int k = 0; k < KS_exp; ++k) 
		{
	  		CipherTextQ *ct = malloc(sizeof *ct);    
	  		ct->b = -old_sk[i]*j*KS_table[k] + Sample(Chi2);
	  		for (int l = 0; l < n; ++l) 
	  		{
	    		ct->a[l] = rand();
	    		ct->b += ct->a[l] * new_sk[l];
	  		}
	  		res[i][j][k] = ct;
		}
}
  
void KeySwitch(CipherTextQ* res, const SwitchingKey K, const CipherTextQN& ct) {
  for (int k = 0; k < n; ++k) 
    res->a[k] = 0;
    
  res->b = ct.b;
  for (int i = 0; i < N; ++i) 
  {
    uZmodQ a = -ct.a[i];
    for (int j = 0; j < KS_exp; ++j, a /= KS_base) 
    {
      uZmodQ a0 = a % KS_base;
      for (int k = 0; k < n; ++k)
        res->a[k] -= (K[i][a0][j])->a[k];
      res->b -= (K[i][a0][j])->b;
    }
  } 
}
  


