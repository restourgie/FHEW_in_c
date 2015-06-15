#include "FFT.h"
#include "params.h"
#include <stdio.h>


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifndef DOUBLE_N
#define DOUBLE_N N*2
#endif

void BitInvert(double complex data[]){
  int i,mv,k,rev;
  double complex temp;

  for(i = 1; i<(DOUBLE_N);i++){//run through all index 1 to N
    k = i;
    mv = (DOUBLE_N)/2;
    rev = 0;
    while(k > 0){//Invert the index
      if ((k % 2) > 0)
              rev = rev + mv;
            k = k / 2;
            mv = mv / 2;
    }
    if(i < rev){
      temp = data[rev];
      data[rev] = data[i];
      data[i] = temp;
    }
  }
}

void CalcFFT(double complex data[], int sign){
  BitInvert(data);
  //variables for the fft
  unsigned long mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

  //Danielson-Lanzcos routine
  mmax=1;
  while ((DOUBLE_N) > mmax) {
    istep=mmax << 1;
    theta=-sign*(2*M_PI/istep);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=0;m<mmax;++m) {
      for (i=m;i<(DOUBLE_N);i+=istep) {
        j=i+mmax;
        tempr = wr * creal(data[j])-wi * cimag(data[j]);
        tempi = wr * cimag(data[j])+wi * creal(data[j]);
        data[j] = data[i] - (tempr + tempi*I);
        data[i] += (tempr + tempi*I);
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
  //end of the algorithm
}

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
void FFTforward(Ring_FFT res, Ring_ModQ val) {
    double complex data[DOUBLE_N];
    for(int k=0;k<N;++k){
      data[k] = val[k] + 0.0*I;
      data[k+N] = 0.0;
    }
    CalcFFT(data,1);

    for(int k=0; k < N2-1; ++k){;
      res[k] = data[2*k+1];
    }
    res[N2-1] = (double complex) 0.0;
}

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
void FFTbackward(Ring_ModQ res, Ring_FFT val){
  double complex data[DOUBLE_N];
  for(int k = 0;k < N2-1; ++k){
    data[2*k+1] = val[k]/N;
    data[2*k] = 0.0;
    data[DOUBLE_N-(2*k+1)] = conj(val[k])/N;
    data[DOUBLE_N-(2*k+2)] = 0.0;
  }
  data[2*N2] = 0.0;

  CalcFFT(data,-1);
  for(int k=0; k < N; ++k)
    res[k] = (long int) round(creal(data[k]));
}


