#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "../mul.h"
#include "negacyclic.h"


double ***wortel;

void print256_num(__m256d var) 
{
    double *v64val = (double*) &var;
    printf("%f %f %f %f\n", v64val[0], v64val[1],v64val[2],v64val[3]);
}

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
{	int m = n/2;
	if(n > 8)
	{
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
	else if(n == 8)
	{	
		level +=2;
	    __m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag,sub_real,sub_imag;
		real_x = _mm256_load_pd(x->real+lo);
		imag_x = _mm256_load_pd(x->imag+lo);
		real_y = _mm256_load_pd(x->real+lo+m);
		imag_y = _mm256_load_pd(x->imag+lo+m);
		//LAYER 2
		//we need to collect all the left and the right parts
		//real part
		sub_real = _mm256_unpacklo_pd(real_x,real_y);
		sub_imag = _mm256_unpackhi_pd(real_x,real_y);

		real_y = _mm256_add_pd(sub_real,sub_imag);
		real_x = _mm256_sub_pd(sub_real,sub_imag);

		//imag part
		sub_real = _mm256_unpacklo_pd(imag_x,imag_y);
		sub_imag = _mm256_unpackhi_pd(imag_x,imag_y);

		imag_y = _mm256_add_pd(sub_real,sub_imag);
		imag_x = _mm256_sub_pd(sub_real,sub_imag);

		//Now do complex inverse computation with roots of unity
		real_twid = _mm256_setr_pd(wortel[0][level][lo/2],wortel[0][level][(lo+m)/2],wortel[0][level][(lo+2)/2],wortel[0][level][(lo+m+2)/2]);
    	imag_twid = _mm256_setr_pd(wortel[2][level][lo/2],wortel[2][level][(lo+m)/2],wortel[2][level][(lo+2)/2],wortel[2][level][(lo+m+2)/2]);
	    temp_real = _mm256_mul_pd(imag_x,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_x,real_twid);

		//TEMP_imag = ad + bc
		imag_x = _mm256_fmadd_pd(real_x,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_x = _mm256_fmsub_pd(real_x,real_twid,temp_real);

		//LAYER 4
		--level;
		//real part
		sub_real = _mm256_permute2f128_pd(real_y,real_x,0x20);
		sub_imag = _mm256_permute2f128_pd(real_y,real_x,0x31);

		real_y = _mm256_add_pd(sub_real,sub_imag);
		real_x = _mm256_sub_pd(sub_real,sub_imag);
		//imag part
		sub_real = _mm256_permute2f128_pd(imag_y,imag_x,0x20);
		sub_imag = _mm256_permute2f128_pd(imag_y,imag_x,0x31);

		imag_y = _mm256_add_pd(sub_real,sub_imag);
		imag_x = _mm256_sub_pd(sub_real,sub_imag);
		//mult
		real_twid = _mm256_setr_pd(wortel[0][level][lo/4],wortel[0][level][(lo+m)/4],wortel[0][level][lo/4],wortel[0][level][(lo+m)/4]);
	    imag_twid = _mm256_setr_pd(wortel[2][level][lo/4],wortel[2][level][(lo+m)/4],wortel[2][level][lo/4],wortel[2][level][(lo+m)/4]);
	    temp_real = _mm256_mul_pd(imag_x,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_x,real_twid);
		//TEMP_imag = ad + bc
		imag_x = _mm256_fmadd_pd(real_x,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_x = _mm256_fmsub_pd(real_x,real_twid,temp_real);
		
		//LAYER 8
		--level;
		//real part
		temp_real = _mm256_unpacklo_pd(real_y,real_x);
		temp_imag = _mm256_unpackhi_pd(real_y,real_x);

		real_x = _mm256_add_pd(temp_real,temp_imag);
		real_y = _mm256_sub_pd(temp_real,temp_imag);
		//imag part
		sub_real  = _mm256_unpacklo_pd(imag_y,imag_x);
		sub_imag  = _mm256_unpackhi_pd(imag_y,imag_x);

		imag_x = _mm256_add_pd(sub_real,sub_imag);
		imag_y = _mm256_sub_pd(sub_real,sub_imag);
		//mult
	    real_twid = _mm256_set1_pd(wortel[0][level][lo/n]);
    	imag_twid = _mm256_set1_pd(wortel[2][level][lo/n]);
	    temp_real = _mm256_mul_pd(imag_y,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_y,real_twid);

		//TEMP_imag = ad + bc
		imag_y = _mm256_fmadd_pd(real_y,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_y = _mm256_fmsub_pd(real_y,real_twid,temp_real);

		real_x =_mm256_permute4x64_pd(real_x,0xd8);
		real_y =_mm256_permute4x64_pd(real_y,0xd8);
		imag_x =_mm256_permute4x64_pd(imag_x,0xd8);
		imag_y =_mm256_permute4x64_pd(imag_y,0xd8);

		_mm256_store_pd(x->real+lo,real_x);
		_mm256_store_pd(x->imag+lo,imag_x);
		_mm256_store_pd(x->real+lo+m,real_y);
		_mm256_store_pd(x->imag+lo+m,imag_y);
	}
}

void recursive_phi(cplx_ptr *x,int n,int lo,int level)
{	
  int m = n/2;	
  if(n > 8)
  {
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
  else if(n == 8)
  {
    __m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag,sub_real,sub_imag;
    real_twid = _mm256_set1_pd(wortel[0][level][lo/n]);
    imag_twid = _mm256_set1_pd(wortel[1][level][lo/n]);

	//(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
	real_x = _mm256_load_pd(x->real+lo);
	imag_x = _mm256_load_pd(x->imag+lo);
	real_y = _mm256_load_pd(x->real+lo+m);
	imag_y = _mm256_load_pd(x->imag+lo+m);
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

	//START LAYER 4
	//We know that real_x = |a|b|c|d| imag_x = |ai|bi|ci|di| real_y = |e|f|g|h| imag_y = |ei|fi|gi|hi| 
	//We need to twidle (c+ci),(d+di) and (g+gi),(h+hi)
	++level;

	sub_real = _mm256_permute2f128_pd(real_x,real_y,0x31);
	sub_imag = _mm256_permute2f128_pd(imag_x,imag_y,0x31);
	real_twid = _mm256_setr_pd(wortel[0][level][lo/4],wortel[0][level][lo/4],wortel[0][level][(lo+m)/4],wortel[0][level][(lo+m)/4]);
    imag_twid = _mm256_setr_pd(wortel[1][level][lo/4],wortel[1][level][lo/4],wortel[1][level][(lo+m)/4],wortel[1][level][(lo+m)/4]);

    temp_real = _mm256_mul_pd(sub_imag,imag_twid);
    temp_imag = _mm256_mul_pd(sub_imag,real_twid);
	//TEMP_real = ac - bd
	temp_real = _mm256_fmsub_pd(sub_real,real_twid,temp_real);
	//TEMP_imag = ad + bc
	temp_imag = _mm256_fmadd_pd(sub_real,imag_twid,temp_imag);
	
	//REAL PART
	//get abef
	sub_real = _mm256_permute2f128_pd(real_x,real_y,0x20);
	//abef  
	//cdgh-
	sub_imag = _mm256_sub_pd(sub_real,temp_real);
	//abef
	//cdgh+
	sub_real = _mm256_add_pd(sub_real,temp_real);

	// //NEEDED TO COMPLETE LAYER 4
	// real_x = _mm256_permute2f128_pd(sub_imag,sub_real,0x02);
	// real_y = _mm256_permute2f128_pd(sub_imag,sub_real,0x13);
	//PREPARE REALS FOR LAYER 2
	//STORE ALL FACTORS THAT NEED TO BE MULT WITH ROOTS OF UNITY IN REAL_X
	real_x = _mm256_unpackhi_pd(sub_real,sub_imag);
	real_y = _mm256_unpacklo_pd(sub_real,sub_imag);


	//IMAG PART
	//get ai bi ei fi
	sub_real = _mm256_permute2f128_pd(imag_x,imag_y,0x20);
	//abef  
	//cdgh-
	sub_imag = _mm256_sub_pd(sub_real,temp_imag);
	//abef
	//cdgh+
	sub_real = _mm256_add_pd(sub_real,temp_imag);

	//NEEDED TO COMPLETE LAYER 4
	// imag_x = _mm256_permute2f128_pd(sub_imag,sub_real,0x02);
	// imag_y = _mm256_permute2f128_pd(sub_imag,sub_real,0x13);
	//PREPARE IMAGS FOR LAYER 2
	//STORE ALL FACTORS THAT NEED TO BE MULT WITH ROOTS OF UNITY IN IMAG_X
	//STORE ALL NON TWIDLE IMAGS IN IMAG_Y 
	imag_x = _mm256_unpackhi_pd(sub_real,sub_imag);
	imag_y = _mm256_unpacklo_pd(sub_real,sub_imag);

	// // //START LAYER 2!!
	++level;
	real_twid = _mm256_setr_pd(wortel[0][level][lo/2],wortel[0][level][(lo+2)/2],wortel[0][level][(lo+m)/2],wortel[0][level][(lo+m+2)/2]);
    imag_twid = _mm256_setr_pd(wortel[1][level][lo/2],wortel[1][level][(lo+2)/2],wortel[1][level][(lo+m)/2],wortel[1][level][(lo+m+2)/2]);

    temp_real = _mm256_mul_pd(imag_x,imag_twid);
    temp_imag = _mm256_mul_pd(imag_x,real_twid);

	//TEMP_real = ac - bd
	temp_real = _mm256_fmsub_pd(real_x,real_twid,temp_real);
	//TEMP_imag = ad + bc
	temp_imag = _mm256_fmadd_pd(real_x,imag_twid,temp_imag);

	sub_real = _mm256_sub_pd(real_y,temp_real);
	sub_imag = _mm256_add_pd(real_y,temp_real); 

	temp_real = _mm256_unpacklo_pd(sub_imag,sub_real);
	sub_real  = _mm256_unpackhi_pd(sub_imag,sub_real);

	real_x = _mm256_permute2f128_pd(temp_real,sub_real,0x20);
	real_y = _mm256_permute2f128_pd(temp_real,sub_real,0x31);

	sub_real = _mm256_sub_pd(imag_y,temp_imag);
	sub_imag = _mm256_add_pd(imag_y,temp_imag);

	temp_real = _mm256_unpacklo_pd(sub_imag,sub_real);
	sub_real  = _mm256_unpackhi_pd(sub_imag,sub_real);

	imag_x = _mm256_permute2f128_pd(temp_real,sub_real,0x20);
	imag_y = _mm256_permute2f128_pd(temp_real,sub_real,0x31);

	_mm256_store_pd(x->real+lo,real_x);
	_mm256_store_pd(x->imag+lo,imag_x);
	_mm256_store_pd(x->real+lo+m,real_y);
	_mm256_store_pd(x->imag+lo+m,imag_y);	
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