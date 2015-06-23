#include <complex.h>
#include <fftw3.h>
#include <iostream>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>


#define N 1024
#define N2 (N/2+1)


typedef int32_t ZmodQ;
typedef ZmodQ Ring_ModQ[N];
typedef fftw_complex Ring_FFT[N2];

using namespace std;

double *in;
fftw_complex *out;
fftw_plan plan_fft_forw, plan_fft_back;


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const int DOUBLE_N = N*2;

/*****************************************************************************************************************
*
*                   Support code
*
******************************************************************************************************************/

//THIS IS NEEDED TO CHECK DOUBLES EQUALITY
int adjust_num(double num) {
    double low_bound = 1e7;
    double high_bound = low_bound*10;
    double adjusted = num;
    int is_negative = (num < 0);
    if(num == 0) {
        return 0;
    }
    if(is_negative) {
        adjusted *= -1;
    }
    while(adjusted < low_bound) {
        adjusted *= 10;
    }
    while(adjusted >= high_bound) {
        adjusted /= 10;
    }
    if(is_negative) {
        adjusted *= -1;
    }
    //define int round(double) to be a function which rounds
    //correctly for your domain application.
    return round(adjusted);
}

inline void wait_on_enter()
{
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
}

/*****************************************************************************************************************
*
*                   My FFT
*
******************************************************************************************************************/

void BitInvert(complex double data[]){
  int i,mv,k,rev;
  complex double temp;

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

void CalcFFT(complex double data[], int sign){
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
void FFTforward_my_version(Ring_FFT res, Ring_ModQ val) {
    complex double data[DOUBLE_N];
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
void FFTbackward_my_version(Ring_ModQ res, Ring_FFT val){
  complex double data[DOUBLE_N];
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

/*****************************************************************************************************************
*
*                   FFTW
*
******************************************************************************************************************/

void FFTforward_fftw_version(Ring_FFT res, const Ring_ModQ val) {
  for (int k = 0; k < N; ++k)	{
    in[k] = (double) (val[k]);
    in[k+N] = 0.0;			
  }
  fftw_execute(plan_fft_forw); 
  for (int k = 0; k < N2; ++k) 
    res[k] = (double complex) out[2*k+1];				
}

void FFTbackward_FFTW_version(Ring_ModQ res, const Ring_FFT val){
  for (int k = 0; k < N2; ++k) {
    out[2*k+1] = (double complex) val[k]/N;
    out[2*k]   = (double complex) 0;
  }
  fftw_execute(plan_fft_back); 
  for (int k = 0; k < N; ++k) 
    res[k] = (long int) round(in[k]);
}


/*****************************************************************************************************************
*
*                   Benchmark functions
*
******************************************************************************************************************/


void FFTforward(Ring_FFT res, const Ring_ModQ val) {
  complex_double result[N2];
  FFTforward_my_version(result,val);
  FFTforward_fftw_version(res,val);

  for(int i=0;i<N2;++i){
    if(adjust_num(creal(res[i]))!= adjust_num(result[i][0])){
      cout << "Alert real doesnt match! Index = " << i << "   Result FFTW Real = " << creal(res[i]) << " Complex = " << cimag(res[i]) << "  Result my FFT Real = " << result[i][0] << " Complex = " << result[i][1] << endl;
      cout << "Int Form of Real. FFTW Real = " << adjust_num(creal(res[i])) << " MY FFT Real = " << adjust_num(result[i][0]) << endl;
    }
    if(adjust_num(cimag(res[i])) != adjust_num(result[i][1])){
      cout << "Alert imaginary doenst match! Index = " << i << "   Result FFTW Real = " << creal(res[i]) << " Complex = " << cimag(res[i]) << "  Result my FFT Real = " << result[i][0] << " Complex = " << result[i][1] << endl;
    }
  }
}

void FFTbackward(Ring_ModQ res,const Ring_FFT val){
  Ring_ModQ result;
  FFTbackward_my_version(result,val);
  FFTbackward_FFTW_version(res,val);
  for(int i =0; i<N;++i){
    if(res[i] != result[i])
      cout << "Index: " << i <<"FFTW Result = " << res[i] << "   My Result = " << result[i] << endl;
  }
}


/*****************************************************************************************************************
*
*                   Main
*
******************************************************************************************************************/

int main(){
  /*TODO: GENERATE A LOT OF FFT's afterwards test the time */
  return 0;
}