#include "LWE.h"
#include "params.h"


void KeyGen(SecretKey sk) {
	KeyGenRestart:
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
  
}

void Encrypt(CipherText* ct, const SecretKey sk, int m) {

}

int Decrypt(const SecretKey sk, const CipherText& ct) {
 
}

void DecryptDetail(const SecretKey sk, const CipherText& ct) {

}

    
int round_qQ(ZmodQ v) {
  
}
  
void ModSwitch(CipherText* ct, const CipherTextQ& c) {

}
  
void SwitchingKeyGen(SwitchingKey res, const SecretKey new_sk, const SecretKeyN old_sk) {

}
  
void KeySwitch(CipherTextQ* res, const SwitchingKey K, const CipherTextQN& ct) {

}
  

