#include "FFT.h"
#include "params.h"


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifndef DOUBLE_N
#define DOUBLE_N N*2
#endif

//typedef double complex_double[2];
void BitInvert(complex_double data[]){
  int i,mv,k,rev;
  complex_double temp;

  for(i = 1; i<(DOUBLE_N);i++){//run through all index 1 to N
    k = i;
    mv = N; //THIS USED TO BE mv = N/2 HOWEVER NOW WE USE TWICE AS BIG FFT SO AND THUS IT BECOMES N*2/2 WHICH MEANS WE CAN JUST PUT N DO NOT FORGET TO CHANGE AFTER
    rev = 0;
    while(k > 0){//Invert the index
      if ((k % 2) > 0)
              rev = rev + mv;
            k = k / 2;
            mv = mv / 2;
    }
    if(i < rev){

      temp[0] = data[rev][0];
      temp[1] = data[rev][1];
      data[rev][0] = data[i][0];
      data[rev][1] = data[i][1];
      data[i][0] = temp[0];
      data[i][1] = temp[1];

    }
  }
}

void CalcFFT(complex_double data[], int sign){
  BitInvert(data);
  //variables for the fft
  unsigned long mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

  //Danielson-Lanzcos routine
  mmax=1;
  while (DOUBLE_N > mmax) {
    istep=mmax << 1;
    theta=sign*(2*M_PI/istep);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=0;m<mmax;++m) {
      for (i=m;i<(DOUBLE_N);i+=istep) {
        j=i+mmax;
        tempr=wr * data[j][0]-wi * data[j][1];
        tempi=wr * data[j][1]+wi * data[j][0];

        data[j][0]= data[i][0]-tempr;
        data[j][1]= data[i][1]-tempi;

        data[i][0] += tempr;
        data[i][1] += tempi;
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
    complex_double data[DOUBLE_N];
    for(int k=0;k<N;++k){
      data[k][0] = val[k];
      data[k][1] = 0.0;
      data[k+N][0] = 0.0;
      data[k+N][1] = 0.0;
    } 
    CalcFFT(data,1);
    for(int k=0; k < N2-1; ++k){
      res[k][0] = data[2*k+1][0];
      res[k][1] = data[2*k+1][1];
    }
    res[N2][0] = 0.0;
    res[N2][1] = 0.0;
}

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
void FFTbackward(Ring_ModQ res, Ring_FFT val){
  complex_double data[DOUBLE_N];
  for(int k = 0;k < N2-1; ++k){
    data[2*k+1][0] = val[k][0]/(double)N; 
    data[2*k+1][1] = val[k][1]/(double)N;
    data[2*k][0] = 0.0;
    data[2*k][1] = 0.0;

    data[DOUBLE_N-(2*k+1)][0] = val[k][0]/(double)N;
    data[DOUBLE_N-(2*k+1)][1] = -(val[k][1]/(double)N);
    data[DOUBLE_N-(2*k+2)][0] = 0.0;
    data[DOUBLE_N-(2*k+2)][1] = 0.0;
  }
  data[2*N2][0] = 0.0;
  data[2*N2][1] = 0.0;

  CalcFFT(data,-1);
  for(int k=0; k < N; ++k)
    res[k] = (long int) round(data[k][0]);
}


