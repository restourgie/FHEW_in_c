#include "fft.h"
#include <math.h>
#include <stdio.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define N 1024
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

void CalcFFT(double complex data[], int sign){
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

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
// void FFTforward(double complex *r, const uint32_t *x)
// {
//     double complex data[DOUBLE_N];
//     int k;
//     for(k=0;k<N;++k){
//       data[k] = x[k] + 0.0*I;
//       data[k+N] = 0.0;
//     }
//     CalcFFT(data,1);

//     for(k=0; k < N2-1; ++k){;
//       r[k] = data[2*k+1];
//     }
//     r[N2-1] = (double complex) 0.0;
// }

// void FFTbackward(uint32_t *r,  const double complex *x)
// {
//   double complex data[DOUBLE_N];
//   int k;
//   for(k = 0;k < N2-1; ++k){
//     data[2*k+1] = x[k]/N;
//     data[2*k] = 0.0;
//     data[DOUBLE_N-(2*k+1)] = conj(x[k])/N;
//     data[DOUBLE_N-(2*k+2)] = 0.0;
//   }
//   data[2*N2] = 0.0;

//   CalcFFT(data,-1);
//   for(k=0; k < N; ++k)
//     r[k] = (long int) round(creal(data[k]));
// }

void FFTsetup() {
  in = (double*) fftw_malloc(sizeof(double) * 2*N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N + 2));
  plan_fft_forw = fftw_plan_dft_r2c_1d(2*N, in, out,  FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_c2r_1d(2*N, out, in,  FFTW_PATIENT);
}

void FFTforward(double complex *r, const uint32_t *x) {
  for (int k = 0; k < N; ++k) {
    in[k] = (double) (x[k]);
    in[k+N] = 0.0;      
  }
  fftw_execute(plan_fft_forw); 
  for (int k = 0; k < N2; ++k) 
    r[k] = (double complex) out[2*k+1];       
}

void FFTbackward(uint32_t *r, const double complex *x){
  for (int k = 0; k < N2; ++k) {
    out[2*k+1] = (double complex) x[k]/N;
    out[2*k]   = (double complex) 0;
  }
  fftw_execute(plan_fft_back); 
  for (int k = 0; k < N; ++k) 
    r[k] = (long int) round(in[k]);
}
