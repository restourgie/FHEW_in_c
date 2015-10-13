#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fft_negacyc_lut.h"
#include "../mul.h"

double complex **wortel;

void init_wortel(int n,int lo,int level,double complex root)
{	
	if(n > 1){
		int m = n/2;
		wortel[level][lo/n] = root;
		++level;
		init_wortel(m,lo,level,csqrt(root));
		init_wortel(m,lo+m,level,csqrt(-root));
	}
}


void init_negacyc()
{
	int loga = log2(CPLXDIM);
	wortel = malloc(loga*sizeof(*wortel));
	int j = 1;
	for (int i = 0; i < loga; ++i)
	{
		wortel[i] = malloc(j*sizeof(double complex));
		j = j <<1;
	}
	init_wortel(CPLXDIM,0,0,csqrt(I));
}

/******************************************************************
*
* SMART COMPLEX MULTIPLICATION
*
******************************************************************/
void inverse_phi_lut(double complex *x,int n,int lo,int level)
{	
	if(n > 1){
		int m = n/2;
		inverse_phi_lut(x,m,lo,level+1);
		inverse_phi_lut(x,m,lo+m,level+1);
		double complex temp;
			for(int i=lo;i<m+lo;++i){
				temp = x[i] + x[i+m];
				x[i+m] =  (x[i] - x[i+m]) * conj(wortel[level][lo/n]);
				x[i] = temp;
			}
	}
}

void recursive_phi_lut(double complex *x,int n,int lo,int level)
{	
  	if(n > 1){
    double complex temp;
    int m = n/2;
    // printf("LEVEL = %d, lo = %d n = %d\n",level,lo,n);    
    // printf("Root: ( %f + I * %f)\n",creal(wortel[level][lo/n]),cimag(wortel[level][lo/n]));
	    for(int i=lo; i < lo+m;++i){
	      temp = wortel[level][lo/n] * x[i+m];
		      //phiprime
	      x[i+m] = x[i] - temp;
		      //phi
	      x[i] = x[i] + temp;
	    }
	++level;
    recursive_phi_lut(x,m,lo,level);
    recursive_phi_lut(x,m,lo + m,level);
  }
}