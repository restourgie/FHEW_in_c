#ifndef LWE_H
#define LWE_H

#include "params.h"
#include "distrib.h"

  
  typedef struct { // n-dimensional LWE ciphertext modulo q
    int a[n]; // array of n integers modulo q 
    int b;    // integer modulo q
  } CipherText;

  typedef struct { // n-dimensional LWE ciphertext modulo Q
    ZmodQ a[n];
    ZmodQ b;
  } CipherTextQ;
  
  typedef struct { // N-dimensional LWE ciphertext modulo Q
    ZmodQ a[N];
    ZmodQ b;
  } CipherTextQN;
  
  typedef int SecretKey[n]; // n-dimensional LWE secret key 
  typedef int SecretKeyN[N]; // N-dimensional LWE secret key 
  typedef CipherTextQ SwitchingKey[N][KS_base][KS_exp]; //SwitchingKey => CipherTextQ[1024][25][7]
  
  SecretKey* LWEKeyGen();
  SecretKeyN* KeyGenN();
  // Generate key material (SwitchingKey) required by KeySwitch to transform 
  // LWE encryptions under old_sk into LWE encryptions under new_sk
  void SwitchingKeyGen(SwitchingKey* res,SecretKey *new_sk, SecretKeyN *old_sk);

  void Encrypt(CipherText* ct, SecretKey* sk, int m);
  int Decrypt(SecretKey sk, CipherText ct);
  
 
  void KeySwitch(CipherTextQ* res, SwitchingKey* Key, CipherTextQN* ct);

  // Changes an LWE ciphertext modulo Q into an LWE ciphertext modulo q
  void ModSwitch(CipherText* ct, CipherTextQ* c); 
  // int round_qQ(ZmodQ v);


// For debbugging purpose you can use the following 
//  void DecryptDetail(const SecretKey sk, const CipherText& ct); 





#endif
