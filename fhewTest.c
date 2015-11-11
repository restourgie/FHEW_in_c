#include "LWE.h"
#include "FHEW.h"
#include "distrib.h"
#include <stdlib.h>

void help(char* cmd) {
  printf("\nusage: %s  n\n\n",cmd); 
  printf("  Generate a secret key sk and evaluation key ek, and repeat the following test n times:\n");
  printf("   - generate random bits b1,b2,b3,b4\n");
  printf("   - compute ciphertexts c1, c2, c3 and c4 encrypting b1, b2, b3 and b4  under sk\n");
  printf("   - homomorphically compute the encrypted (c1 NAND c2) NAND (c3 NAND c4) \n");
  printf("   - decrypt all the intermediate results and check correctness \n");
  printf("\n If any of the tests fails, print ERROR and stop immediately.\n\n");
  exit(0);
}

int main(int argc, char *argv[]) {
  if (argc != 2) help(argv[0]);
  int count = atoi(argv[1]); 

  printf("Setting up FHEW \n");
  Setup();
  printf("Generating secret key ... ");
  SecretKey LWEsk;
  LWEKeyGen(LWEsk);
  printf(" Done.\n");
  printf("Generating evaluation key ... this may take a while ... ");
  EvalKey EK;

  EK.BSkey = (BootstrappingKey*) malloc(sizeof(BootstrappingKey));
  if(EK.BSkey == NULL) {
    fprintf(stderr, "BAD BAD BAD!\n");
    return -1;
  }
  
  EK.KSkey = (SwitchingKey*) malloc(sizeof(SwitchingKey));
  if(EK.KSkey == NULL) {
    fprintf(stderr, "EVEN WORSE!\n");
    return -1;
  }

  FHEWKeyGen(&EK, LWEsk);
  printf(" Done.\n\n");
  printf("Testing homomorphic NAND %d times.\n",count); 
  printf("Circuit shape : (a NAND b) NAND (c NAND d)\n\n");



  for (int i = 0; i < count; ++i) {
    int a,b,c,d,ab,cd,abcd;
    CipherText a_e,b_e,c_e,d_e,ab_e,cd_e,abcd_e;

    a = rand()%2;  
    b = rand()%2;
    c = rand()%2;
    d = rand()%2;
    Encrypt(&a_e, LWEsk, a);
    Encrypt(&b_e, LWEsk, b);
    Encrypt(&c_e, LWEsk, c);
    Encrypt(&d_e, LWEsk, d);
  
    int testa,testb,testc,testd;
    testa = Decrypt(LWEsk,&a_e);
    testb = Decrypt(LWEsk,&b_e);
    testc = Decrypt(LWEsk,&c_e);
    testd = Decrypt(LWEsk,&d_e);
    if(a != testa)
      printf("ERROR DECRYPT FAULT at A\n");
    if(b != testb)
      printf("ERROR DECRYPT FAULT at B\n");
    if(c != testc)
      printf("ERROR DECRYPT FAULT at C\n");
    if(d != testd)
      printf("ERROR DECRYPT FAULT at D\n");

    printf("Enc(%d)  NAND  Enc(%d)  =  ",a,b);
    HomNAND(&ab_e, &EK, &a_e, &b_e);
    ab = Decrypt(LWEsk, &ab_e);
    printf("Enc(%d)\n",ab);

    if(1 - a*b != ab){ 
      printf("ERROR at ab iteration %d\n",i); 
      exit(1); 
    }

    printf("Enc(%d)  NAND  Enc(%d)  =  ",c,d);
    HomNAND(&cd_e, &EK, &c_e, &d_e);
    cd = Decrypt(LWEsk, &cd_e);
    printf("Enc(%d)\n",cd);

    if(1 - c*d != cd){ 
      printf("ERROR at cd iteration %d\n",i); 
      exit(1); 
    }

    printf("Enc(%d)  NAND  Enc(%d)  =  ",ab,cd);
    HomNAND(&abcd_e, &EK, &ab_e, &cd_e);
    abcd = Decrypt(LWEsk, &abcd_e);
    printf("Enc(%d)\n",abcd);

    if(1 - ab*cd != abcd){ 
      printf("ERROR at abcd iteration %d\n",i); 
      exit(1); 
    }
    printf("\nNext iteration\n");
}

  printf("\nPassed all tests!\n\n");
}


