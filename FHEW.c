#include "FHEW.h"
#include "FFT.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

// int compare(const void * elem1, const void * elem2)
// { 
//   uint64_t x = *((uint64_t*)elem1);
//   uint64_t y = *((uint64_t*)elem2);
//   if(x > y) 
//     return 1;
//   else if(x < y) 
//     return -1;
//   return 0;
// }

// uint64_t rdtsc()
// {
//     unsigned int lo,hi;
//     __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
//     return ((uint64_t)hi << 32) | lo;
// }

typedef Ring_ModQ ct_ModQ[K2][2];   // Ciphertext in coefficient form. ct_ModQ => Ring_ModQ[6][2] => Ring_ModQ[6][2] =>  ZmodQ[6][2][1024] => int32_t[6][2][1024]
typedef Ring_ModQ dct_ModQ[K2][K2]; // Decomposed Ciphertext in coeff form. dct_ModQ => Ring_ModQ[6][6] => ZmodQ[6][6][1024] => int32_t[6][6][1024]
typedef Ring_FFT  dct_FFT[K2][K2];  // Decomposed Ciphertext in FFT form. dct_FFT => Ring_FFT[6][6] => complex double[6][6][513]

Ring_FFT t_TestMSB;

void Setup(){
  FFTsetup();
	Ring_ModQ tmsb;
  tmsb[0]=-1;
	for (int i = 1; i < N; ++i)
		tmsb[i]=1;
	FFTforward(t_TestMSB,tmsb);
}

/*************************************************************************
*                                                                        *
*                           ENCRYPTION                                   *
*                                                                        *
*************************************************************************/

void FHEWencrypt(ct_FFT ct, Ring_FFT sk_FFT, int m) {
    Ring_FFT ai;
    ct_ModQ res;
    int mm = (((m % q) + q) % q) * (2*N/q);// Reduce mod q (dealing with negative number as well)
    int sign = 1;
    
    if (mm >= N) { 
      mm -= N; 
      sign = -1; 
    }
    
    for (int i = 0; i < K2; ++i) 
    {
      for (int k = 0; k < N; ++k) 
        res[i][0][k] = random_int(); // % Q //
      
      FFTforward(ai, res[i][0]);
      
      for (int k = 0; k < N2; ++k)
        ai[k] = ai[k] * sk_FFT[k];

      FFTbackward(res[i][1], ai);
      
      for (int k = 0; k < N; ++k) 
        res[i][1][k] += Sample_1(Chi1); // Add error [a,as+e]
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

/*************************************************************************
*                                                                        *
*                           KEY GENERATION                               *
*                                                                        *
*************************************************************************/
//BootstrappingKey[500][23][2] => ct_FFT[500][23][2] => ct_FFT[500][23][2][6][2] => complex_double[500][23][2][6][2][513] => double[500][23][2][6][2][513][2]
void FHEWKeyGen(EvalKey* EK, SecretKey LWEsk){
  SecretKeyN FHEWsk;
  KeyGenN(FHEWsk);
  SwitchingKeyGen((*EK->KSkey),LWEsk, FHEWsk); 

  Ring_FFT FHEWskFFT;
  FFTforward(FHEWskFFT,FHEWsk);

  printf("\nStarting GENERATION of BSkey\n");
  for (int i = 0; i < n; ++i){
    printf("i is now : %d\n",i);
    for (int j = 1; j < BS_base; ++j)
      for (int k = 0; k < BS_exp; ++k) 
        FHEWencrypt( ((*EK->BSkey)[i][j][k]), FHEWskFFT, LWEsk[i] * j * BS_table[k] ); 
  }
  printf("finished GENERATION of BSkey\n");
}

void AddToACC(ct_FFT ACC, ct_FFT C) {
	ct_ModQ ct;		//int32_t[6][2][1024]
  dct_ModQ dct;	//int32_t[6][6][1024]
  dct_FFT dctFFT;	//Complex double[6][6][513]
  for (int i = 0; i < K2; ++i)
    for (int j = 0; j < 2; ++j)
	     FFTbackward(ct[i][j], ACC[i][j]);
  for (int i = 0; i < K2; ++i)
    for (int j = 0; j < 2; ++j)
	     for (int k = 0; k < N; ++k) 
      	{
        		ZmodQ t = ct[i][j][k] * v_inverse;
        		for (int l = 0; l < K; ++l) 
        		{
          		ZmodQ r = (t << g_bits_32[l]) >> g_bits_32[l];
      		    t = (t-r) >> g_bits[l];
      		    dct[i][j+2*l][k] = r;
        		}
      	}
  // uint64_t start,end; 
  // uint64_t cycles[1000];    
  // for (int count = 0; count < 1000; ++count)
  // { 
  //   // printf("count = %d\n",count );
  //   start = rdtsc();  

  for (int i = 0; i < K2; ++i)
    for (int j = 0; j < K2; ++j)
       FFTforward(dctFFT[i][j], dct[i][j]);

  //   end = rdtsc();
  //   cycles[count] = end - start;
  // }
  // qsort(cycles,sizeof(cycles)/sizeof(*cycles),sizeof(*cycles),compare);
  // printf("Median: %llu\n",cycles[499]);

  for (int i = 0; i < K2; ++i)     
    for (int j = 0; j < 2; ++j)
	     for (int k = 0; k < N2; ++k) 
      	{    
            ACC[i][j][k] = 0.0 + 0.0 * I;
        		for (int l = 0; l < K2; ++l)
          		ACC[i][j][k] += dctFFT[i][l][k] * C[l][j][k];//HERE WE ADD *ACC TO THE COMPLEX MULTIPLICATION OF dctFFT[i][l][k] and C[l][j][k]
        
      	}

}

void InitializeACC(ct_FFT ACC, int m) {
	ct_ModQ res;										//int32_t[6][2][1024]
  int mm = (((m % q) + q) % q) * (2*N/q);             // Reduce mod q (dealing with negative number as well)
  int sign = 1;
  if(mm >= N) 
  { 
  	mm -= N; 
  	sign = -1; 
  }
  
  for (int i = 0; i < K2; ++i)
    for (int j = 0; j < 2; ++j)
	    for (int k = 0; k < N; ++k)
  		  res[i][j][k]=0;
  
  for (int i = 0; i < K; ++i) 
  {
    res[2*i  ][0][mm] += sign*vgprime[i]; // Add G Multiple
    res[2*i+1][1][mm] += sign*vgprime[i]; // [a,as+e] + X^m *G
  }
  
  for (int i = 0; i < K2; ++i)
    for (int j = 0; j < 2; ++j)
	    FFTforward(ACC[i][j], res[i][j]);
}


CipherTextQN* MemberTest(Ring_FFT t, ct_FFT C) {
    Ring_FFT temp;
    Ring_ModQ temp_ModQ;
    
    CipherTextQN* ct = malloc(sizeof(CipherTextQN));
    for (int i = 0; i < N2; ++i)
	  	temp[i] = conj(C[1][0][i] * t[i]);  // Compute t*a //temp will be THE COMPLEX MULTIPLICATION OF C[i][l][k] and t[l][j][k]
  
    FFTbackward(temp_ModQ, temp);

    for (int i = 0; i < N; ++i) 
      ct->a[i] = temp_ModQ[i];
  	
    for (int i = 0; i < N2; ++i){
  		temp[i] = C[1][1][i] * t[i];//temp will be THE COMPLEX MULTIPLICATION OF C[i][l][k] and t[l][j][k]
      // printf("C[1][1][%d] = %f + i %f\n",i,creal(C[1][1][i]),cimag(C[1][1][i]));
    }
  	
    FFTbackward(temp_ModQ, temp);
    ct->b = v+temp_ModQ[0];	
    return ct; 
}

void HomNAND(CipherText* res, EvalKey* EK, CipherText* ct1, CipherText* ct2) 
{
	  CipherText e12;
    
    for (int i = 0; i < n; ++i)
      e12.a[i] = (2*q - (ct1->a[i] + ct2->a[i])) % q;
    
    e12.b  =  (13 * q / 8) - (ct1->b + ct2->b) % q;

    ct_FFT ACC;
    InitializeACC(ACC, (e12.b + q/4) % q);
    
    for (int i = 0; i < n; ++i) 
    {
      int a = (q - e12.a[i] % q) % q;
      for (int k = 0; k < BS_exp; ++k, 	a /= BS_base) 
      {
		    int a0 = a % BS_base;
		    if(a0) 
			   AddToACC(ACC, ((*EK->BSkey)[i][a0][k]));	
      }
    }
    CipherTextQN* eQN = MemberTest(t_TestMSB,ACC);
    CipherTextQ *eQ = malloc(sizeof(CipherTextQ));
    KeySwitch(eQ, *(EK->KSkey), eQN);
    free(eQN);
    ModSwitch(res, eQ);
    free(eQ);
}

