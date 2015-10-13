#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "fft_negacyc.h"

/******************************************************************
*
* SMART COMPLEX MULTIPLICATION
*
******************************************************************/
void inverse_phi(double complex *x,int n,int lo,double complex root)
{	
	if(n > 1){
		int m = n/2;
		inverse_phi(x,m,lo,csqrt(root));
		inverse_phi(x,m,lo+m,csqrt(-root));
		double complex temp;
			for(int i=lo;i<m+lo;++i){
				temp = x[i] + x[i+m];
				x[i+m] =  (x[i] - x[i+m]) * conj(root);
				x[i] = temp;
			}
	}
}

void recursive_phi(double complex *x,int n,int lo,double complex root)
{	
  	if(n > 1){
    double complex temp;
    int m = n/2;
    // printf("n = %d, lo = %d\n",n,lo );
    // printf("Root: ( %f + I * %f)\n",creal(root),cimag(root));
	    for(int i=lo; i < lo+m;++i){
	      temp = root * x[i+m];
		      //phiprime
	      x[i+m] = x[i] - temp;
		      //phi
	      x[i] = x[i] + temp;
	    }
    recursive_phi(x,m,lo,csqrt(root));
    recursive_phi(x,m,lo + m,csqrt(-root));
  }
}