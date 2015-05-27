#include "FFT.h"
#include "params.h"
#include <stdio.h>


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifndef DOUBLE_N
#define DOUBLE_N N*2
#endif

void BitInvert(complex_double data[]){
  int i,mv,k,rev;
  complex_double temp;

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
    //printf("\n\n****************Starting FFT****************\n\n");
    complex_double data[DOUBLE_N];
	double t;
    for(int k=0;k<N;++k){
      // printf("Index %d value: %d\n",k,val[k]);
      data[k][0] = val[k];
      data[k][1] = 0.0;
      data[k+N][0] = 0.0;
      data[k+N][1] = 0.0;
    } 
    CalcFFT(data,1);
    // printf("\n\n****************TOTAL RESULT FFT****************\n\n");
    // for(int i =0; i< DOUBLE_N; ++i){
    //   printf("Index: %d real: %f Imag: %f\n",i,data[i][0],data[i][1]);
    // }
    for(int k=0; k < N2-1; ++k){
//	printf("here: %d\n",k);
	t = data[2*k+1][0];
        res[k][0] = t;
	t = data[2*k+1][1];
      	res[k][1] = t;
    }
    res[N2][0] = 0.0;
    res[N2][1] = 0.0;
    // printf("\n\n****************Result FFT****************\n\n");
    // for(int i=0; i < N2; ++i){
    //   printf("Index %d Real: %f Imag: %f\n",i,res[i][0],res[i][1]);
    // }
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


