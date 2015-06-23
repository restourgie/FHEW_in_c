#include <stdio.h>
#include <complex.h>
#include <math.h>


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define N 8
#ifndef DOUBLE_N
#define DOUBLE_N 8
#endif

void wait_on_enter()
{
    printf("Press ENTER to continue\n");
    getchar();
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
void FFTforward(double complex res[], int val[]) {
  double complex data[N];
  
  //data[0] =-1 + I*1;
  for(int i=0;i<N/2;++i)
    data[i] = 1 + I*1;
  for(int i=N/2;i<N;++i)
    data[i] = 0 + I*0;
  
  CalcFFT(data,1);

  double complex temp;
  for(int i=0;i<N-1;++i){
    printf("index i: %d\n",i);
    temp = (cos(2*M_PI*i/N) + I*sin(2*M_PI*i/N));
    printf("temp is: (%f,%f)\n",creal(temp),cimag(temp));
    printf("data[i] = (%f,%f) data[N-i] = (%f,%f)\n",creal(data[i]),cimag(data[i]),creal(data[N-i]),cimag(data[N-i]));
    res[i] = 0.5*(data[i] + conj(data[N-i])) - I*0,5*(data[i] - conj(data[N-i]))*temp;
    printf("the result is: (%f,%f)\n",creal(res[i]),cimag(res[i]));
    //wait_on_enter();
  }
  res[7] = (data[7] + conj(data[0]))/2 + I*((data[7] - conj(data[0]))/2);

}

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
void FFTbackward(double complex data[], int val[]){


  CalcFFT(data,-1);
  for(int k=0; k < N; ++k)
    val[k] = creal(data[k]);
}

int main()
{

  double complex data[N];
  int buf[1];

  FFTforward(data,buf);

  for(int i=0;i<N;++i){
    printf("index i:%d (%f,%f)\n",i,creal(data[i]),cimag(data[i]));
  }
}


