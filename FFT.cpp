#include <complex.h>
#include <fftw3.h>
#include "FFT.h"
#include "params.h"
#include <iostream>


using namespace std;

double *in;
fftw_complex *out;
fftw_plan plan_fft_forw, plan_fft_back;


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const int DOUBLE_N = N*2;


inline void wait_on_enter()
{
    std::string dummy;
    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
}

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
    // cout << "istep " << istep << " theta " << theta << " wtemp " << wtemp << " wpr " << wpr << " wpi " << wpi << " wr " << wr << " wi " << wi << endl;
    for (m=0;m<mmax;++m) {
      for (i=m;i<(DOUBLE_N);i+=istep) {
        j=i+mmax;
        tempr = wr * creal(data[j])-wi * cimag(data[j]);
        tempi = wr * cimag(data[j])+wi * creal(data[j]);
        // cout << "tempr: " << tempr << " tempi: " << tempi << endl;
        // cout << "index j " << j << " creal(data[j]) " << creal(data[j]) << " cimag(data[j]) " << cimag(data[j]) << endl;
        // cout << "index i " << i << " creal(data[i]) " << creal(data[i]) << " cimag(data[i]) " << cimag(data[i]) << endl;
        // cout << "\nafter\n"; 
        data[j] = data[i] - (tempr + tempi*I);
        data[i] += (tempr + tempi*I);
        // cout << "index j " << j << " creal(data[j]) " << creal(data[j]) << " cimag(data[j]) " << cimag(data[j]) << endl;
        // cout << "index i " << i << " creal(data[i]) " << creal(data[i]) << " cimag(data[i]) " << cimag(data[i]) << endl;
        // wait_on_enter();
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      // cout << "wr: " << wr << " wi: " << wi << " wtemp: " << wtemp << endl;
      // wait_on_enter();
    }
    mmax=istep;
  }
  //end of the algorithm
}

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
void FFTforward_my_version(complex double res[],const Ring_ModQ val) {
    complex double data[DOUBLE_N];
    for(int k=0;k<N;++k){
      // cout << "Index " << k << " value: " << val[k] << endl;
      data[k] = val[k] + 0.0*I;
      data[k+N] = 0.0;
    }
    CalcFFT(data,1);
    // for(int i=0;i<(N*2);++i){
    //   cout << "Index " << i << " Real: " << creal(data[i]) << " Imag: " << cimag(data[i]) << endl;
    //   wait_on_enter();
    // }

    for(int k=0; k < N2-1; ++k){
      //  printf("here: %d\n",k);
      res[k] = data[2*k+1];
    }
    res[N2-1] = (double complex) 0.0;

    // for(int i=0; i <= N; ++i){
    //   cout << "FrontIndex: " << i << " Real: " << creal(data[i]) << " Imag: " << cimag(data[i]);
    //   if(i != 0)
    //     cout << "  BackIndex: " << (DOUBLE_N-i) << " Real: " << creal(data[(DOUBLE_N-i)]) << " Imag: " << cimag(data[(DOUBLE_N-i)]);
    //   cout << endl;
    // }
    // wait_on_enter();
}
  
void FFTsetup() {
  in = (double*) fftw_malloc(sizeof(double) * 2*N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N + 2));
  plan_fft_forw = fftw_plan_dft_r2c_1d(2*N, in, out,  FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_c2r_1d(2*N, out, in,  FFTW_PATIENT);
}
  
void FFTforward(Ring_FFT res, const Ring_ModQ val) {
  for (int k = 0; k < N; ++k) {
    in[k] = (double) (val[k]);
    in[k+N] = 0.0;      
  }
  fftw_execute(plan_fft_forw); 
  for (int k = 0; k < N2; ++k) 
    res[k] = (double complex) out[2*k+1];       
}

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

// void FFTforward(Ring_FFT res, const Ring_ModQ val) {
//   complex double result[N2];
//   FFTforward_my_version(result,val);
//   FFTforward_fftw_version(res,val);
//   cout << "\n\n************************************************* RESULTS OF FFT Forward*************************************************\n\n";
//   for(int i=0;i<N2;++i){
//     cout << "index i: " << i << " FFT_REAL : my " << creal(result[i]) << " ,fftw " << creal(res[i]) << " FFT_IMAG : my" << cimag(result[i]) << " , fftw" << cimag(res[i]) << endl;
//     wait_on_enter();
//     // if(adjust_num(creal(res[i]))!= adjust_num(result[i][0])){
//     //   cout << "Alert real doesnt match! Index = " << i << "   Result FFTW Real = " << creal(res[i]) << " Complex = " << cimag(res[i]) << "  Result my FFT Real = " << result[i][0] << " Complex = " << result[i][1] << endl;
//     //   cout << "Int Form of Real. FFTW Real = " << adjust_num(creal(res[i])) << " MY FFT Real = " << adjust_num(result[i][0]) << endl;
//     //   wait_on_enter();
//     // }
//     // if(adjust_num(cimag(res[i])) != adjust_num(result[i][1])){
//     //   cout << "Alert imaginary doenst match! Index = " << i << "   Result FFTW Real = " << creal(res[i]) << " Complex = " << cimag(res[i]) << "  Result my FFT Real = " << result[i][0] << " Complex = " << result[i][1] << endl;
//     //   wait_on_enter();
//     // }
//   }
//   wait_on_enter();
// }


void FFTbackward_FFTW_version(Ring_ModQ res, const Ring_FFT val){
  cout << "\n\n*************************************************STARTING FFTW Backward*************************************************\n\n";
  for (int k = 0; k < N2; ++k) {
    out[2*k+1] = (double complex) val[k]/N;
    cout << "Filling position: " << (2*k+1) << " With real value: " << creal(out[2*k+1]) << " and complex value: " << cimag(out[2*k+1]) << endl;
    out[2*k]   = (double complex) 0;
  }
  fftw_execute(plan_fft_back); 
  for (int k = 0; k < N; ++k) 
    res[k] = (long int) round(in[k]);
}

// // Ring_FFT => complex_double[513] => double[513][2]
// // Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
// void FFTbackward_my_version(Ring_ModQ res,const Ring_FFT val){
//   cout << "\n\n*************************************************STARTING MY FFT Backward*************************************************\n\n";
//   complex_double data[DOUBLE_N];
//   for(int k = 0;k < N2-1; ++k){
//     cout << "Forwards Filling position: " << (2*k+1) << " With real value: " << creal(val[k])/(double)N << " and complex value: " << (cimag(val[k])/(double)N) << endl;
//     data[2*k+1][0] = creal(val[k])/(double)N; 
//     data[2*k+1][1] = cimag(val[k])/(double)N;
//     data[2*k][0] = 0.0;
//     data[2*k][1] = 0.0;

//     data[DOUBLE_N-(2*k+1)][0] = creal(val[k])/(double)N;
//     data[DOUBLE_N-(2*k+1)][1] = -(cimag(val[k]))/(double)N;
//     data[DOUBLE_N-(2*k+2)][0] = 0.0;
//     data[DOUBLE_N-(2*k+2)][1] = 0.0;
//   }
//   data[2*N2][0] = 0.0;
//   data[2*N2][1] = 0.0;


//   CalcFFT(data,-1);
//   for(int k=0; k < N; ++k)
//     res[k] = (long int) round(data[k][0]);
// }

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
void FFTbackward_my_version(Ring_ModQ res, const Ring_FFT val){
  cout << "\n\n*************************************************STARTING MY FFT Backward*************************************************\n\n";
  complex double data[DOUBLE_N];
  for(int k = 0;k < N2-1; ++k){
    data[2*k+1] = val[k]/N;
    data[2*k] = 0.0;
    data[DOUBLE_N-(2*k+1)] = conj(val[k])/N;
    data[DOUBLE_N-(2*k+2)] = 0.0;
  }
  data[2*N2] = 0.0;

  // for(int i=0; i <= N; ++i){
  //     cout << "F_Index: " << i; 
  //     if(i != 0)
  //       cout << "\x1b[32mBackIndex :" << (DOUBLE_N-i) << "\x1b[0m"; 
  //     cout << " Real: " << data[i][0] ;
  //     if(i != 0)
  //       cout << " \x1b[32mReal: " << data[(DOUBLE_N-i)][0] << "\x1b[0m";
  //     cout << "  Imag: " << data[i][1];
  //     if(i != 0)
  //       cout << " \x1b[32mImag: " << data[(DOUBLE_N-i)][1] << "\x1b[0m";
  //     cout << endl;
  // }
  // wait_on_enter();

  CalcFFT(data,-1);
  for(int k=0; k < N; ++k)
    res[k] = (long int) round(creal(data[k]));
    // res[k] = (long int) round(data[k][0]);
}

void FFTbackward(Ring_ModQ res,const Ring_FFT val){
  Ring_ModQ result;
  FFTbackward_my_version(result,val);
  FFTbackward_FFTW_version(res,val);
  cout << "\n\n************************************************* RESULTS OF FFT Backward*************************************************\n\n";
  for(int i =0; i<N;++i){
    // if(res[i] != result[i])
      cout << "Index: " << i <<"FFTW Result = " << res[i] << "   My Result = " << result[i] << endl;
      wait_on_enter();
  }
  wait_on_enter();
}
