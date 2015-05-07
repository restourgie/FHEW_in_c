#ifndef PARAM_H
#define PARAM_H

//#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdint.h>

const int n = 500;

const int N = 1024;
const volatile int N2 = N/2+1;

const int K = 3;
const int K2 = 6;

const volatile long int Q = (long int) 1 << 32;
// #define Q 1 << 32
int q = 512;
int q2 = 256;

typedef int32_t ZmodQ;
typedef uint32_t uZmodQ;

ZmodQ v = (1 << 29) +1;     // Q/8 +1
ZmodQ v_inverse = 3758096385; // 1/v mod Q

ZmodQ vgprime[3] = {v, v<<11, v<<22};
int g_bits[3] = {11, 11, 10};
int g_bits_32[3] = {21, 21, 22};


int KS_base = 25;
int KS_exp = 7;
ZmodQ KS_table[7] = {1,
			   25,
			   25*25,
			   25*25*25,
			   25*25*25*25,
			   25*25*25*25*25,
			   25*25*25*25*25*25};

int BS_base = 23;
int BS_exp = 2;
int BS_table[2] = {1,23};

typedef ZmodQ Ring_ModQ[N]; //Ring_ModQ => ZmodQ[1024] => int32_t[1024]
typedef double complex_double[2];
typedef complex_double Ring_FFT[N2]; //Ring_FFT => complex_double[513] => double[513][2]

#endif
