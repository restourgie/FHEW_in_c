#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "../mul.h"
#include "negacyclic.h"


double ***wortel;

void init_wortel(int n,int lo,int level,double complex root)
{	
	if(n > 1){
		int m = n/2;
		wortel[0][level][lo/n] = creal(root);
		wortel[1][level][lo/n] = cimag(root);
		wortel[2][level][lo/n] = -cimag(root);
		++level;
		init_wortel(m,lo,level,csqrt(root));
		init_wortel(m,lo+m,level,csqrt(-root));
	}
}


void init_negacyc()
{	
	wortel = malloc(3*sizeof(**wortel));
	int loga = log2(CPLXDIM);
	wortel[0] = malloc(loga*sizeof(*wortel));
	wortel[1] = malloc(loga*sizeof(*wortel));
	wortel[2] = malloc(loga*sizeof(*wortel));
	int j = 1;
	for (int i = 0; i < loga; ++i)
	{
		wortel[0][i] = malloc(j*sizeof(double));
		wortel[1][i] = malloc(j*sizeof(double));
		wortel[2][i] = malloc(j*sizeof(double));
		j = j <<1;
	}
	init_wortel(CPLXDIM,0,0,csqrt(I));
}

/******************************************************************
*
* SMART COMPLEX MULTIPLICATION
*
******************************************************************/
void inverse_phi(cplx_ptr *x,int n,int lo,int level)
{	
	if(n > 4)
	{
		int m = n/2;
		inverse_phi(x,m,lo,level+1);
		inverse_phi(x,m,lo+m,level+1);
	    __m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag;
    	real_twid = _mm256_set1_pd(wortel[0][level][lo/n]);
    	imag_twid = _mm256_set1_pd(wortel[2][level][lo/n]);
		for(int i=lo;i<m+lo;i+=4)
		{	
			real_x = _mm256_load_pd(x->real+i);
		    imag_x = _mm256_load_pd(x->imag+i);
		    real_y = _mm256_load_pd(x->real+i+m);
		    imag_y = _mm256_load_pd(x->imag+i+m);

		    temp_real = real_x;
		    temp_imag = imag_x;

		    real_x = _mm256_add_pd(temp_real,real_y);
	  		imag_x = _mm256_add_pd(temp_imag,imag_y);

	  		real_y = _mm256_sub_pd(temp_real,real_y);
	  		imag_y = _mm256_sub_pd(temp_imag,imag_y);

	  		//TEMP_real = bd
		    temp_real = _mm256_mul_pd(imag_y,imag_twid);
		    //TEMP_imag = bc
		    temp_imag = _mm256_mul_pd(imag_y,real_twid);

		    //imag_y = ad + bc
			imag_y = _mm256_fmadd_pd(real_y,imag_twid,temp_imag);
		    //real_y = ac - bd
		    real_y = _mm256_fmsub_pd(real_y,real_twid,temp_real);
			

			_mm256_store_pd(x->real+i,real_x);
		    _mm256_store_pd(x->imag+i,imag_x);
		    _mm256_store_pd(x->real+i+m,real_y);
		    _mm256_store_pd(x->imag+i+m,imag_y);
		}
	}
	else if(n > 1)
	{
		int m = n/2;
		inverse_phi(x,m,lo,level+1);
		inverse_phi(x,m,lo+m,level+1);
		double temp_real,temp_imag;
		double a,b,c,d;
		for(int i=lo;i<m+lo;++i)
		{
			temp_real = x->real[i];
			temp_imag = x->imag[i];

			x->real[i] = temp_real + x->real[i+m];
			x->imag[i] = temp_imag + x->imag[i+m]; 
			
			temp_real = temp_real - x->real[i+m];
			temp_imag = temp_imag - x->imag[i+m]; 

			//(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
			a = temp_real;
		    b = temp_imag;
		    c = wortel[0][level][lo/n];
		    d = wortel[2][level][lo/n];

		    x->real[i+m] = (a*c) - (b*d);
      		x->imag[i+m] = (a*d) + (b*c);
		}
	}
}

void recursive_phi(cplx_ptr *x,int n,int lo,int level)
{	
  if(n > 4)
  {
    int m = n/2;
    __m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag;
    real_twid = _mm256_set1_pd(wortel[0][level][lo/n]);
    imag_twid = _mm256_set1_pd(wortel[1][level][lo/n]);
    for(int i=lo; i < lo+m;i+=4)
    {
	  //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
	  real_x = _mm256_load_pd(x->real+i);
      imag_x = _mm256_load_pd(x->imag+i);
      real_y = _mm256_load_pd(x->real+i+m);
      imag_y = _mm256_load_pd(x->imag+i+m);
      //TEMP_real = bd
      temp_real = _mm256_mul_pd(imag_y,imag_twid);
      //TEMP_imag = bc
      temp_imag = _mm256_mul_pd(imag_y,real_twid);

      //TEMP_real = ac - bd
      temp_real = _mm256_fmsub_pd(real_y,real_twid,temp_real);
	  //TEMP_imag = ad + bc
	  temp_imag = _mm256_fmadd_pd(real_y,imag_twid,temp_imag);

	  real_y = _mm256_sub_pd(real_x,temp_real);
	  imag_y = _mm256_sub_pd(imag_x,temp_imag);

	  real_x = _mm256_add_pd(real_x,temp_real);
	  imag_x = _mm256_add_pd(imag_x,temp_imag);

	  _mm256_store_pd(x->real+i,real_x);
      _mm256_store_pd(x->imag+i,imag_x);
      _mm256_store_pd(x->real+i+m,real_y);
      _mm256_store_pd(x->imag+i+m,imag_y);
	}
	++level;
    recursive_phi(x,m,lo,level);
    recursive_phi(x,m,lo + m,level);
  }
  else if(n > 1)
  {
    double temp_real,temp_imag;
    int m = n/2;
    double a,b,c,d;
    // printf("LEVEL = %d, lo = %d n = %d\n",level,lo,n);    
    // printf("Root: ( %f + I * %f)\n",creal(wortel[level][lo/n]),cimag(wortel[level][lo/n]));
	    for(int i=lo; i < lo+m;++i){
	      //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
	      a = x->real[i+m];
	      b = x->imag[i+m];
	      c = wortel[0][level][lo/n];
	      d = wortel[1][level][lo/n];

	      temp_real = (a*c) - (b*d);
	      temp_imag = (a*d) + (b*c);
		      //phiprime
	      x->real[i+m] = x->real[i] - temp_real;
	      x->imag[i+m] = x->imag[i] - temp_imag;
		      //phi
	      x->real[i] = x->real[i] + temp_real;
	      x->imag[i] = x->imag[i] + temp_imag;
	    }
	++level;
    recursive_phi(x,m,lo,level);
    recursive_phi(x,m,lo + m,level);
  }
}

void phi_forward(cplx_ptr *x,const ring_t *ring)
{
  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    x->real[i] = ring->v[i];
    x->imag[i] = ring->v[j];
    ++j;
  }
  recursive_phi(x, CPLXDIM,0,0);
}

void phi_backward(cplx_ptr *x, ring_t *ring)
{
  inverse_phi(x, CPLXDIM,0,0);

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    ring->v[i] = x->real[i];
    ring->v[j] = x->imag[i];
    ++j;
  }
}