#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "normal_fft.h"
#include "../mul.h"

//TAKES AS INPUT THE ARRAY, THE SIZE, THE CURRENT LOCATION AND THE ROOT
void recursive_FFT(double complex *x,int n,int lo,double complex root)
{ 

  if(n > 1){
    double complex temp;
    int m = n/2;
    // printf("\nN is at the moment: %d\n",n);
    // printf("Root: ( %f + I * %f)\n",creal(root),cimag(root));
    for(int i=lo; i < lo+m;++i){
      temp = root * x[i+m];
      // printf("lo = %d, i = %d, m = %d\n",lo,i,m);
      // printf("temp: ( %f + I * %f)\n",creal(temp),cimag(temp) );
        //phiprime
      x[i+m] = x[i] - temp;
        //phi
      x[i] = x[i] + temp;
    }
    // print_complex(x,REALDIM);
    recursive_FFT(x,m,lo,csqrt(root));
    if(root == 1){
      root = I;
      recursive_FFT(x,m,lo+m,root);
    }
    else  
      recursive_FFT(x,m,lo + m,csqrt(-root));
  }
}

void inverse_FFT(double complex *x,int n,int lo,double complex root)
{ 
  if(n > 1){
    int m = n/2;
    inverse_FFT(x,m,lo,csqrt(root));
    if(root == 1){
      inverse_FFT(x,m,lo+m,I);
    }
    else
      inverse_FFT(x,m,lo+m,csqrt(-root));

    // printf("\nN is at the moment: %d\n",n);
    // printf("lo = %d, m = %d\n",lo,m);
    // printf("Root: ( %f + I * %f)\n",creal(root),cimag(root));
    double complex temp;
    for(int i=lo;i<m+lo;++i){
      // printf("x[i]= ( %f + I * %f) x[i+m] = ( %f + I * %f)\n",creal(x[i]),cimag(x[i]),creal(x[i+m]),cimag(x[i+m]) );
      temp = x[i] + x[i+m];
      x[i+m] =  ((x[i] - x[i+m]) * conj(root));

      x[i] = temp;
      // printf("x[i]= ( %f + I * %f) x[i+m] = ( %f + I * %f)\n",creal(x[i]),cimag(x[i]),creal(x[i+m]),cimag(x[i+m]) );
    }
    //print_complex(x,REALDIM);
  }
}