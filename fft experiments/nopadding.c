#include <complex.h>
#include <fftw3.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>


#define N 1024
#define N2 (N/2+1)
#define NR_TESTS 100


typedef int32_t ZmodQ;
typedef ZmodQ Ring_ModQ[N];
typedef fftw_complex Ring_FFT[N2];

double *in;
fftw_complex *out;
fftw_plan plan_fft_forw, plan_fft_back;


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const int DOUBLE_N = N;


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

int comp (const void * elem1, const void * elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

inline void wait_on_enter()
{
    printf("Press ENTER to continue\n");
    getchar();
}

/*****************************************************************************************************************
*
*                   Rosetta code
*
******************************************************************************************************************/

typedef double complex cplx;
 
void _fft(cplx buf[], cplx out[], int n, int step)
{
  if (step < n) {
    _fft(out, buf, n, step * 2);
    _fft(out + step, buf + step, n, step * 2);
 
    for (int i = 0; i < n; i += 2 * step) {
      cplx t = cexp(-I * M_PI * i / n) * out[i + step];
      buf[i / 2]     = out[i] + t;
      buf[(i + n)/2] = out[i] - t;
    }
  }
}
 
void fft(cplx buf[], int n)
{
  cplx out[n];
  for (int i = 0; i < n; i++) out[i] = buf[i];
 
  _fft(buf, out, n, 1);
}

void show(const char * s, cplx buf[]) {
  printf("%s", s);
  for (int i = 0; i < N; i++)
    if (!cimag(buf[i]))
      printf("%g ", creal(buf[i]));
    else
      printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
}

/*****************************************************************************************************************
*
*                   My FFT
*
******************************************************************************************************************/
void my_setup(){

}


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
void FFTforward_my_version(Ring_FFT res, Ring_ModQ val) {
    double complex data[DOUBLE_N];
    for(int k=0;k<N;++k){
      data[k] = val[k] + 0.0*I;
      // data[k+N] = 0.0;
    }
    CalcFFT(data,1);
    printf("\n\n MY FFT result\n\n");
    for(int i=0;i<N;++i){
      printf("Index: %d Real: %f, Imag: %f\n",i,creal(data[i]),cimag(data[i]));
    }
    printf("\n\n****** end *****\n\n");

    for(int k=0; k < N2-1; ++k){;
      res[k] = data[2*k+1];
    }
    res[N2-1] = (double complex) 0.0;
}

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
// void FFTbackward_my_version(Ring_ModQ res, Ring_FFT val){
//   double complex data[DOUBLE_N];
//   for(int k = 0;k < N2-1; ++k){
//     data[2*k+1] = val[k]/N;
//     data[2*k] = 0.0;
//     data[DOUBLE_N-(2*k+1)] = conj(val[k])/N;
//     data[DOUBLE_N-(2*k+2)] = 0.0;
//   }
//   data[2*N2] = 0.0;

//   CalcFFT(data,-1);
//   for(int k=0; k < N; ++k)
//     res[k] = (long int) round(creal(data[k]));
// }

/*****************************************************************************************************************
*
*                   FFTW
*
******************************************************************************************************************/

void FFTW_setup(){
  in = (double*) fftw_malloc(sizeof(double) * 2*N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N + 2));
  plan_fft_forw = fftw_plan_dft_r2c_1d(2*N, in, out,  FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_c2r_1d(2*N, out, in,  FFTW_PATIENT);
}

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

static __inline__ unsigned long GetCC(void)
{
  unsigned a, d; 
  asm volatile("rdtsc" : "=a" (a), "=d" (d)); 
  return ((unsigned long)a) | (((unsigned long)d) << 32); 
}


void FFTtest() {
  double complex result[N2];
  //Ring_FFT res[NR_TESTS];
  //Ring_ModQ val[NR_TESTS];
  Ring_FFT res;
  // unsigned long long t[NR_TESTS];

   Ring_ModQ tmsb;   
   tmsb[0]=-1;    
   for (int i = 1; i < N; ++i)    
    tmsb[i]=1; 

  cplx buf[N];
  buf[0] = -1;
  for(int i =1;i<N;++i)
    buf[i] = 1;
  show("Data: ",buf);
  fft(buf,N);
  show("\nFFT: ",buf);

  FFTforward_my_version(result,tmsb);
  FFTforward_fftw_version(res,tmsb);

  for(int i=0;i<N2;++i){
   printf("Index = %d FFT Real = %f, %f Complex = %f, %f\n",i,creal(res[i]),creal(result[i]),cimag(res[i]),cimag(result[i]));
    // if(adjust_num(creal(res[i]))!= adjust_num(cimag(result[i]))){
    //   printf("Alert real doesnt match! Index = %d Result FFTW Real = %f Complex = %f Result my FFT Real = %f Complex = %f\n",i,creal(res[i]),cimag(res[i]),creal(result[i]),cimag(result[i]));
    //   printf("Int Form of Real. FFTW Real = %d MY FFT Real = %d\n",adjust_num(creal(res[i])),adjust_num(creal(result[i])));
    // }
    // if(adjust_num(cimag(res[i])) != adjust_num(cimag(result[i]))){
    //   printf("Alert imaginary doenst match! Index = %d Result FFTW Real = %f Complex = %f Result my FFT Real = %f Complex = %f\n",i,creal(res[i]),cimag(res[i]),creal(result[i]),cimag(result[i]));
    // }
  }
}

void FFTbackward(Ring_ModQ res, Ring_FFT val){
  Ring_ModQ result;
  // FFTbackward_my_version(result,val);
  FFTbackward_FFTW_version(res,val);
  for(int i =0; i<N;++i){
    if(res[i] != result[i])
      printf("Index: %d FFTW Result = %d My Result = %d\n",i,res[i],result[i]);
  }
}


/*****************************************************************************************************************
*
*                   Main
*
******************************************************************************************************************/

int main(){
  FFTW_setup();

  FFTtest();

  /*TODO: GENERATE A LOT OF FFT's afterwards test the time */
  return 0;
}