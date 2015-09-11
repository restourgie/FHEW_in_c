#include "fft.h"
#include <math.h>
#include <stdio.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define N 4
#define N2 (N/2+1)

#ifndef DOUBLE_N
#define DOUBLE_N N*2
#endif

double *in;
fftw_complex *out;
fftw_plan plan_fft_forw, plan_fft_back;

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

void CalcFFT(double complex data[], int sign)
{
  BitInvert(data);
  //variables for the fft
  unsigned long mmax,m,j,istep,i;
  double wtemp,theta;
  double complex twiddle,wp,temp;
  //Danielson-Lanczos routine
  mmax=1;
  while ((DOUBLE_N) > mmax) {
    // printf("poof\n");
    istep=mmax << 1;
    theta=-sign*(2*M_PI/istep);
    wtemp=sin(0.5*theta);
    wp = -2.0*wtemp*wtemp + I*sin(theta);
    twiddle = 1.0 + I*0;
    for (m=0;m<mmax;++m) 
    {
      // printf("looping\n");
      for (i=m;i<(DOUBLE_N);i+=istep) 
      {
        // printf("i:%d\n",i);
        j=i+mmax;
        temp = twiddle * data[j];
        data[j] = data[i] - temp;
        data[i] += temp;
      }
      twiddle += twiddle*wp;
    }
    mmax=istep;
  }
  //end of the algorithm
}

// void FFTforward(double complex *r, const uint32_t *x)
// {
//   double complex data[DOUBLE_N];
//   int k;
//   for(k=0;k<N/2;++k){
//     data[k] = x[k] + x[(N/2)+k]*I;
//   }
//   CalcFFT(data,1);

//   for(k=0; k < N2-1; ++k){;
//     r[k] = data[k];
//   }
//   r[N2-1] = (double complex) 0.0;
// }

// ORIGINAL VERSION
void FFTforward(double complex *r, const uint32_t *x)
{
  double complex data[DOUBLE_N];
  int k;

  for(k=0;k<N;++k){
    data[k] = x[k] + 0.0*I;
    data[k+N] = 0.0;
  }
  CalcFFT(data,1);

  for(k=0; k < N2-1; ++k){;
    r[k] = data[2*k+1];
  }
  r[N2-1] = (double complex) 0.0;
}

void CPLX_FFTforward(double complex *r)
{
  CalcFFT(r,1);
}


// //TESTING OTHER TYPE
// void FFTforward(double complex *r, const uint32_t *x)
// {
//   //double complex data[DOUBLE_N];
//   int k;
//   for(k=0;k<N;++k){
//     r[k] = x[k] + 0.0*I;
//     r[k+N] = 0.0;
//   }
//   CalcFFT(r,1);


  
// }

// //TESTING OTHER TYPE
// void FFTbackward(uint32_t *r, double complex *x)
// {
//   CalcFFT(x,-1);
//   int k;
//   for(k=0; k < N; ++k)
//     r[k] = (long int) round(creal(x[k]) - creal(x[k+N]))/(DOUBLE_N);
// }



// void FFTbackward(uint32_t *r,  const double complex *x)
// {
//   double complex data[DOUBLE_N];
//   int k;
//   for(k = 0;k < N2-1; ++k){
//     data[k] = x[k]/N;
//   }
//   CalcFFT(data,-1);
//   for(k=0; k < N/2; ++k){
//     r[k] = (long int) round(creal(data[k]));
//     r[(N/2)+k] = (long int) round(cimag(data[k]));
//   }
// }


// ORIGINAL VERSION
void FFTbackward(uint32_t *r,  const double complex *x)
{
  double complex data[DOUBLE_N];
  int k;
  for(k = 0;k < N2-1; ++k){
    data[2*k+1] = x[k]/N;
    data[2*k] = 0.0;
    data[DOUBLE_N-(2*k+1)] = conj(x[k])/N;
    data[DOUBLE_N-(2*k+2)] = 0.0;
  }
  data[2*N2] = 0.0;

  CalcFFT(data,-1);
  for(k=0; k < N; ++k)
    r[k] = (long int) round(creal(data[k]));
}

void CPLX_FFTbackward(double complex *r)
{
    CalcFFT(r,-1);
}
