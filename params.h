#ifndef PARAM_H
#define PARAM_H

#include <math.h>
#include <stdint.h>

#define n 500

#define N 1024
#define N2 N/2+1

#define K 3
#define K2 6

const long long Q = (long long) 1 << 32;
const int q = 512;
const int q2 = 256;

typedef int32_t ZmodQ;
typedef uint32_t uZmodQ;
const ZmodQ v = (1 << 29) +1;       // Q/8 +1
const ZmodQ v_inverse = 3758096385; // 1/v mod Q

const ZmodQ vgprime[3] = {(long long)(1 << 29) +1, (long long)1<<11, (long long)1<<22};
const int g_bits[3] = {11, 11, 10};
const int g_bits_32[3] = {21, 21, 22};


#define KS_base 25
#define KS_exp 7
const ZmodQ KS_table[7] = {1,
			   25,
			   25*25,
			   25*25*25,
			   25*25*25*25,
			   25*25*25*25*25,
			   25*25*25*25*25*25};

#define BS_base 23
#define BS_exp 2
const int BS_table[2] = {1,23};

typedef ZmodQ Ring_ModQ[N]; //Ring_ModQ => ZmodQ[1024] => int32_t[1024]
typedef double complex_double[2];
typedef complex_double Ring_FFT[N2]; //Ring_FFT => complex_double[513] => double[513][2]

#endif