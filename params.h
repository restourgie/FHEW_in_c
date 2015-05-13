#ifndef PARAM_H
#define PARAM_H

#include <math.h>
#include <stdint.h>

#define n 500

#define N 1024
#define N2 N/2+1

#define K 3
#define K2 6

// const long long Q = (long long) 1 << 32;
#define Q 4294967296LL
#define q 512
#define q2 256

typedef int32_t ZmodQ;
typedef uint32_t uZmodQ;

#define v 536870913LL       // Q/8 +1
#define v_inverse 3758096385LL // 1/v mod Q

extern const ZmodQ vgprime[3];
extern const int g_bits[3];
extern const int g_bits_32[3];


#define KS_base 25
#define KS_exp 7
extern const ZmodQ KS_table[7];

#define BS_base 23
#define BS_exp 2
extern const int BS_table[2];

typedef ZmodQ Ring_ModQ[N]; //Ring_ModQ => ZmodQ[1024] => int32_t[1024]
typedef double complex_double[2];
typedef complex_double Ring_FFT[N2]; //Ring_FFT => complex_double[513] => double[513][2]

#endif