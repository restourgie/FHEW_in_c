#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "../mul.h"
#include "negacyclic.h"


double ***wortel;

void print_complx(const cplx_ptr *x,int N){
  for (int i = 0; i < N; ++i)
    printf("cplxpoly[%d] = %f + i * %f\n",i,x->real[i],x->imag[i]);
  printf("\n");
}

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
		--level;
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
	int j = CPLXDIM/2;
	for (int i = 0; i < loga; ++i)
	{
		posix_memalign((void**)&wortel[0][i],32,j*sizeof(double));
		posix_memalign((void**)&wortel[1][i],32,j*sizeof(double));
		posix_memalign((void**)&wortel[2][i],32,j*sizeof(double));
		j = j>>1;
	}
	init_wortel(CPLXDIM,0,(loga-1),csqrt(I));
}

/******************************************************************
*
* SMART COMPLEX MULTIPLICATION
*
******************************************************************/
void inverse_phi(cplx_ptr *x)
{	
	__m256d temp_real,temp_imag,sub_real,sub_imag,v0_r,v128_r,v256_r
		,v384_r,v0_i,v128_i,v256_i,v384_i,real_twid_L1,
  		imag_twid_L1,real_twid_L2_1,imag_twid_L2_1,real_twid_L2_2,
  		imag_twid_L2_2;
  	int lo;

	for(lo = 0;lo < CPLXDIM;lo +=8)
  	{
		v0_r = _mm256_load_pd(x->real+lo);
		v0_i = _mm256_load_pd(x->imag+lo);
		v256_r = _mm256_load_pd(x->real+lo+4);
		v256_i = _mm256_load_pd(x->imag+lo+4);
		//LAYER 2
		//we need to collect all the left and the right parts
		//real part
		sub_real = _mm256_unpacklo_pd(v0_r,v256_r);
		sub_imag = _mm256_unpackhi_pd(v0_r,v256_r);

		v256_r = _mm256_add_pd(sub_real,sub_imag);
		v0_r = _mm256_sub_pd(sub_real,sub_imag);

		//imag part
		sub_real = _mm256_unpacklo_pd(v0_i,v256_i);
		sub_imag = _mm256_unpackhi_pd(v0_i,v256_i);

		v256_i = _mm256_add_pd(sub_real,sub_imag);
		v0_i = _mm256_sub_pd(sub_real,sub_imag);

		//Now do complex inverse computation with roots of unity
		real_twid_L1 = _mm256_setr_pd(wortel[0][0][lo/2],wortel[0][0][(lo+4)/2],wortel[0][0][(lo+2)/2],wortel[0][0][(lo+6)/2]);
    	imag_twid_L1 = _mm256_setr_pd(wortel[2][0][lo/2],wortel[2][0][(lo+4)/2],wortel[2][0][(lo+2)/2],wortel[2][0][(lo+6)/2]);
	    temp_real = _mm256_mul_pd(v0_i,imag_twid_L1);
    	temp_imag = _mm256_mul_pd(v0_i,real_twid_L1);

		//TEMP_imag = ad + bc
		v0_i = _mm256_fmadd_pd(v0_r,imag_twid_L1,temp_imag);
		//TEMP_real = ac - bd
		v0_r = _mm256_fmsub_pd(v0_r,real_twid_L1,temp_real);

		//LAYER 4
		//real part
		sub_real = _mm256_permute2f128_pd(v256_r,v0_r,0x20);
		sub_imag = _mm256_permute2f128_pd(v256_r,v0_r,0x31);

		v256_r = _mm256_add_pd(sub_real,sub_imag);
		v0_r = _mm256_sub_pd(sub_real,sub_imag);
		//imag part
		sub_real = _mm256_permute2f128_pd(v256_i,v0_i,0x20);
		sub_imag = _mm256_permute2f128_pd(v256_i,v0_i,0x31);

		v256_i = _mm256_add_pd(sub_real,sub_imag);
		v0_i = _mm256_sub_pd(sub_real,sub_imag);
		//mult
		real_twid_L1 = _mm256_setr_pd(wortel[0][1][lo/4],wortel[0][1][(lo+4)/4],wortel[0][1][lo/4],wortel[0][1][(lo+4)/4]);
	    imag_twid_L1 = _mm256_setr_pd(wortel[2][1][lo/4],wortel[2][1][(lo+4)/4],wortel[2][1][lo/4],wortel[2][1][(lo+4)/4]);
	    temp_real = _mm256_mul_pd(v0_i,imag_twid_L1);
    	temp_imag = _mm256_mul_pd(v0_i,real_twid_L1);
		//TEMP_imag = ad + bc
		v0_i = _mm256_fmadd_pd(v0_r,imag_twid_L1,temp_imag);
		//TEMP_real = ac - bd
		v0_r = _mm256_fmsub_pd(v0_r,real_twid_L1,temp_real);
		
		//LAYER 8
		//real part
		temp_real = _mm256_unpacklo_pd(v256_r,v0_r);
		temp_imag = _mm256_unpackhi_pd(v256_r,v0_r);

		v0_r = _mm256_add_pd(temp_real,temp_imag);
		v256_r = _mm256_sub_pd(temp_real,temp_imag);
		//imag part
		sub_real  = _mm256_unpacklo_pd(v256_i,v0_i);
		sub_imag  = _mm256_unpackhi_pd(v256_i,v0_i);

		v0_i = _mm256_add_pd(sub_real,sub_imag);
		v256_i = _mm256_sub_pd(sub_real,sub_imag);
		//mult
	    real_twid_L1 = _mm256_set1_pd(wortel[0][2][lo/8]);
    	imag_twid_L1 = _mm256_set1_pd(wortel[2][2][lo/8]);
	    temp_real = _mm256_mul_pd(v256_i,imag_twid_L1);
    	temp_imag = _mm256_mul_pd(v256_i,real_twid_L1);

		//TEMP_imag = ad + bc
		v256_i = _mm256_fmadd_pd(v256_r,imag_twid_L1,temp_imag);
		//TEMP_real = ac - bd
		v256_r = _mm256_fmsub_pd(v256_r,real_twid_L1,temp_real);

		v0_r =_mm256_permute4x64_pd(v0_r,0xd8);
		v256_r =_mm256_permute4x64_pd(v256_r,0xd8);
		v0_i =_mm256_permute4x64_pd(v0_i,0xd8);
		v256_i =_mm256_permute4x64_pd(v256_i,0xd8);

		_mm256_store_pd(x->real+lo,v0_r);
		_mm256_store_pd(x->imag+lo,v0_i);
		_mm256_store_pd(x->real+lo+4,v256_r);
		_mm256_store_pd(x->imag+lo+4,v256_i);
  	}

  //Merge 32-16
  for(lo = 0;lo < CPLXDIM;lo +=32)
  {	
	//init twiddle for 32
	real_twid_L1 = _mm256_set1_pd(wortel[0][4][lo/32]);
	imag_twid_L1 = _mm256_set1_pd(wortel[2][4][lo/32]);
	//init twiddle for 16
	real_twid_L2_1 = _mm256_set1_pd(wortel[0][3][lo/16]);
	imag_twid_L2_1 = _mm256_set1_pd(wortel[2][3][lo/16]);
	real_twid_L2_2 = _mm256_set1_pd(wortel[0][3][(lo+16)/16]);
	imag_twid_L2_2 = _mm256_set1_pd(wortel[2][3][(lo+16)/16]);

  	for(int offset = 0; offset < 8; offset+=4)
  	{
  		//LOAD A LOT of STUF
	  	//FIRST ROUND WE NEED 0..3,32..35 LEFT SIDE (Real and Imag)
	  	//64..67,96..99 (Real and Imag)
	  	//REAL PART
	  	v0_r   = _mm256_load_pd(x->real+lo+offset);
	  	v128_r = _mm256_load_pd(x->real+lo+offset+8);
	  	v256_r = _mm256_load_pd(x->real+lo+offset+16);
	  	v384_r = _mm256_load_pd(x->real+lo+offset+24);
	  	//IMAG PART
		v0_i  = _mm256_load_pd(x->imag+lo+offset);
	  	v128_i  = _mm256_load_pd(x->imag+lo+offset+8);
	  	v256_i = _mm256_load_pd(x->imag+lo+offset+16);
	  	v384_i = _mm256_load_pd(x->imag+lo+offset+24);

		//WE START WITH LAYER 64
	  	//TWIDDLE 0 - 32
		temp_real = _mm256_sub_pd(v0_r,v128_r);
		temp_imag = _mm256_sub_pd(v0_i,v128_i);

		v0_r = _mm256_add_pd(v0_r,v128_r);
		v0_i = _mm256_add_pd(v0_i,v128_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid_L2_1);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid_L2_1);

		v128_i = _mm256_fmadd_pd(temp_real,imag_twid_L2_1,temp_imag);
		v128_r = _mm256_fmsub_pd(temp_real,real_twid_L2_1,sub_real);
		//TWIDDLE 64 - 96
		temp_real = _mm256_sub_pd(v256_r,v384_r);
		temp_imag = _mm256_sub_pd(v256_i,v384_i);

		v256_r = _mm256_add_pd(v256_r,v384_r);
		v256_i = _mm256_add_pd(v256_i,v384_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid_L2_2);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid_L2_2);

		v384_i = _mm256_fmadd_pd(temp_real,imag_twid_L2_2,temp_imag);
		v384_r = _mm256_fmsub_pd(temp_real,real_twid_L2_2,sub_real);

		//NOW WE start Layer 128
	  	//TWIDDLE 0 - 64
		temp_real = _mm256_sub_pd(v0_r,v256_r);
		temp_imag = _mm256_sub_pd(v0_i,v256_i);

		v0_r = _mm256_add_pd(v0_r,v256_r);
		v0_i = _mm256_add_pd(v0_i,v256_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid_L1);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid_L1);

		v256_i = _mm256_fmadd_pd(temp_real,imag_twid_L1,temp_imag);
		v256_r = _mm256_fmsub_pd(temp_real,real_twid_L1,sub_real);
	  	//TWIDDLE 32 - 96
		temp_real = _mm256_sub_pd(v128_r,v384_r);
		temp_imag = _mm256_sub_pd(v128_i,v384_i);

		v128_r = _mm256_add_pd(v128_r,v384_r);
		v128_i = _mm256_add_pd(v128_i,v384_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid_L1);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid_L1);

		v384_i = _mm256_fmadd_pd(temp_real,imag_twid_L1,temp_imag);
		v384_r = _mm256_fmsub_pd(temp_real,real_twid_L1,sub_real);		


		//STORE ALL RESULTS
		_mm256_store_pd(x->real+lo+offset,v0_r);
		_mm256_store_pd(x->real+lo+offset+8,v128_r);
		_mm256_store_pd(x->real+lo+offset+16,v256_r);
		_mm256_store_pd(x->real+lo+offset+24,v384_r);

		_mm256_store_pd(x->imag+lo+offset,v0_i);
		_mm256_store_pd(x->imag+lo+offset+8,v128_i);
		_mm256_store_pd(x->imag+lo+offset+16,v256_i);
		_mm256_store_pd(x->imag+lo+offset+24,v384_i);
  	}
  }


  //Merge 128-64
  for(lo = 0;lo < CPLXDIM;lo +=128)
  {	
	//init twiddle for 128
	real_twid_L1 = _mm256_set1_pd(wortel[0][6][lo/128]);
	imag_twid_L1 = _mm256_set1_pd(wortel[2][6][lo/128]);
	//init twiddle for 64
	real_twid_L2_1 = _mm256_set1_pd(wortel[0][5][lo/64]);
	imag_twid_L2_1 = _mm256_set1_pd(wortel[2][5][lo/64]);
	real_twid_L2_2 = _mm256_set1_pd(wortel[0][5][(lo+64)/64]);
	imag_twid_L2_2 = _mm256_set1_pd(wortel[2][5][(lo+64)/64]);

  	for(int offset = 0; offset < 32; offset+=4)
  	{	
  		//LOAD A LOT of STUF
	  	//FIRST ROUND WE NEED 0..3,32..35 LEFT SIDE (Real and Imag)
	  	//64..67,96..99 (Real and Imag)
	  	//REAL PART
	  	v0_r   = _mm256_load_pd(x->real+lo+offset);
	  	v128_r = _mm256_load_pd(x->real+lo+offset+32);
	  	v256_r = _mm256_load_pd(x->real+lo+offset+64);
	  	v384_r = _mm256_load_pd(x->real+lo+offset+96);
	  	//IMAG PART
		v0_i  = _mm256_load_pd(x->imag+lo+offset);
	  	v128_i  = _mm256_load_pd(x->imag+lo+offset+32);
	  	v256_i = _mm256_load_pd(x->imag+lo+offset+64);
	  	v384_i = _mm256_load_pd(x->imag+lo+offset+96);

	  	//WE START WITH LAYER 64
	  	//TWIDDLE 0 - 32
		temp_real = _mm256_sub_pd(v0_r,v128_r);
		temp_imag = _mm256_sub_pd(v0_i,v128_i);

		v0_r = _mm256_add_pd(v0_r,v128_r);
		v0_i = _mm256_add_pd(v0_i,v128_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid_L2_1);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid_L2_1);

		v128_i = _mm256_fmadd_pd(temp_real,imag_twid_L2_1,temp_imag);
		v128_r = _mm256_fmsub_pd(temp_real,real_twid_L2_1,sub_real);
		//TWIDDLE 64 - 96
		temp_real = _mm256_sub_pd(v256_r,v384_r);
		temp_imag = _mm256_sub_pd(v256_i,v384_i);

		v256_r = _mm256_add_pd(v256_r,v384_r);
		v256_i = _mm256_add_pd(v256_i,v384_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid_L2_2);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid_L2_2);

		v384_i = _mm256_fmadd_pd(temp_real,imag_twid_L2_2,temp_imag);
		v384_r = _mm256_fmsub_pd(temp_real,real_twid_L2_2,sub_real);

		//NOW WE start Layer 128
	  	//TWIDDLE 0 - 64
		temp_real = _mm256_sub_pd(v0_r,v256_r);
		temp_imag = _mm256_sub_pd(v0_i,v256_i);

		v0_r = _mm256_add_pd(v0_r,v256_r);
		v0_i = _mm256_add_pd(v0_i,v256_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid_L1);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid_L1);

		v256_i = _mm256_fmadd_pd(temp_real,imag_twid_L1,temp_imag);
		v256_r = _mm256_fmsub_pd(temp_real,real_twid_L1,sub_real);
	  	//TWIDDLE 32 - 96
		temp_real = _mm256_sub_pd(v128_r,v384_r);
		temp_imag = _mm256_sub_pd(v128_i,v384_i);

		v128_r = _mm256_add_pd(v128_r,v384_r);
		v128_i = _mm256_add_pd(v128_i,v384_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid_L1);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid_L1);

		v384_i = _mm256_fmadd_pd(temp_real,imag_twid_L1,temp_imag);
		v384_r = _mm256_fmsub_pd(temp_real,real_twid_L1,sub_real);		

		//STORE ALL RESULTS
		_mm256_store_pd(x->real+lo+offset,v0_r);
		_mm256_store_pd(x->real+lo+offset+32,v128_r);
		_mm256_store_pd(x->real+lo+offset+64,v256_r);
		_mm256_store_pd(x->real+lo+offset+96,v384_r);

		_mm256_store_pd(x->imag+lo+offset,v0_i);
		_mm256_store_pd(x->imag+lo+offset+32,v128_i);
		_mm256_store_pd(x->imag+lo+offset+64,v256_i);
		_mm256_store_pd(x->imag+lo+offset+96,v384_i);
  	}
  }



  //init twiddle for 512
  real_twid_L1 = _mm256_set1_pd(wortel[0][8][0]);
  imag_twid_L1 = _mm256_set1_pd(wortel[2][8][0]);
  //init twiddle for 256
  real_twid_L2_1 = _mm256_set1_pd(wortel[0][7][0]);
  imag_twid_L2_1 = _mm256_set1_pd(wortel[2][7][0]);
  real_twid_L2_2 = _mm256_set1_pd(wortel[0][7][1]);
  imag_twid_L2_2 = _mm256_set1_pd(wortel[2][7][1]);

  //FUSING 512-256
  for(int offset = 0; offset < 128; offset+=4)
  {	
	//LOAD A LOT of STUF
  	//FIRST ROUND WE NEED 0..3,128..131 LEFT SIDE (Real and Imag)
  	//256..259,384..387 (Real and Imag)
  	//REAL PART
  	v0_r  = _mm256_load_pd(x->real+offset);
  	v128_r  = _mm256_load_pd(x->real+offset+128);
  	v256_r = _mm256_load_pd(x->real+offset+256);
  	v384_r = _mm256_load_pd(x->real+offset+384);
  	//IMAG PART
	v0_i  = _mm256_load_pd(x->imag+offset);
  	v128_i  = _mm256_load_pd(x->imag+offset+128);
  	v256_i = _mm256_load_pd(x->imag+offset+256);
  	v384_i = _mm256_load_pd(x->imag+offset+384);

  	//WE START WITH LAYER 256
  	//TWIDDLE 0 - 128
	temp_real = _mm256_sub_pd(v0_r,v128_r);
	temp_imag = _mm256_sub_pd(v0_i,v128_i);

	v0_r = _mm256_add_pd(v0_r,v128_r);
	v0_i = _mm256_add_pd(v0_i,v128_i);

	sub_real = _mm256_mul_pd(temp_imag,imag_twid_L2_1);
	temp_imag = _mm256_mul_pd(temp_imag,real_twid_L2_1);

	v128_i = _mm256_fmadd_pd(temp_real,imag_twid_L2_1,temp_imag);
	v128_r = _mm256_fmsub_pd(temp_real,real_twid_L2_1,sub_real);
	//TWIDDLE 256 - 284
	temp_real = _mm256_sub_pd(v256_r,v384_r);
	temp_imag = _mm256_sub_pd(v256_i,v384_i);

	v256_r = _mm256_add_pd(v256_r,v384_r);
	v256_i = _mm256_add_pd(v256_i,v384_i);

	sub_real = _mm256_mul_pd(temp_imag,imag_twid_L2_2);
	temp_imag = _mm256_mul_pd(temp_imag,real_twid_L2_2);

	v384_i = _mm256_fmadd_pd(temp_real,imag_twid_L2_2,temp_imag);
	v384_r = _mm256_fmsub_pd(temp_real,real_twid_L2_2,sub_real);

	//NOW WE start Layer 512
  	//TWIDDLE 0 - 256
	temp_real = _mm256_sub_pd(v0_r,v256_r);
	temp_imag = _mm256_sub_pd(v0_i,v256_i);

	v0_r = _mm256_add_pd(v0_r,v256_r);
	v0_i = _mm256_add_pd(v0_i,v256_i);

	sub_real  = _mm256_mul_pd(temp_imag,real_twid_L1);
	temp_imag = _mm256_mul_pd(temp_real,real_twid_L1);

	v256_i = _mm256_sub_pd(sub_real,temp_imag);
	v256_r = _mm256_add_pd(sub_real,temp_imag);
  	//TWIDDLE 128 - 384
	temp_real = _mm256_sub_pd(v128_r,v384_r);
	temp_imag = _mm256_sub_pd(v128_i,v384_i);

	v128_r = _mm256_add_pd(v128_r,v384_r);
	v128_i = _mm256_add_pd(v128_i,v384_i);

	sub_real  = _mm256_mul_pd(temp_imag,real_twid_L1);
	temp_imag = _mm256_mul_pd(temp_real,real_twid_L1);

	v384_i = _mm256_sub_pd(sub_real,temp_imag);
	v384_r = _mm256_add_pd(sub_real,temp_imag);
	
	//STORE ALL RESULTS
	_mm256_store_pd(x->real+offset,v0_r);
	_mm256_store_pd(x->real+offset+128,v128_r);
	_mm256_store_pd(x->real+offset+256,v256_r);
	_mm256_store_pd(x->real+offset+384,v384_r);

	_mm256_store_pd(x->imag+offset,v0_i);
	_mm256_store_pd(x->imag+offset+128,v128_i);
	_mm256_store_pd(x->imag+offset+256,v256_i);
	_mm256_store_pd(x->imag+offset+384,v384_i);
  }

}

void iterative_phi(cplx_ptr *x)
{
  __m256d temp_real,temp_imag,sub_real,sub_imag,real_twid_L1,
  		imag_twid_L1,real_twid_L2_1,imag_twid_L2_1,real_twid_L2_2,
  		imag_twid_L2_2,v0_r,v128_r,v256_r,v384_r,v0_i,v128_i,v256_i,v384_i;
  int lo;

  //init twiddle for 512
  //ONLY init 1 because REAL = IMAG
  // real_twid_L1 = _mm256_set1_pd(wortel[0][8][0]);
  imag_twid_L1 = _mm256_set1_pd(wortel[1][8][0]);
  //init twiddle for 256
  real_twid_L2_1 = _mm256_set1_pd(wortel[0][7][0]);
  imag_twid_L2_1 = _mm256_set1_pd(wortel[1][7][0]);
  real_twid_L2_2 = _mm256_set1_pd(wortel[0][7][1]);
  imag_twid_L2_2 = _mm256_set1_pd(wortel[1][7][1]);

  //FUSING 512-256
  for(int offset = 0; offset < 128; offset+=4)
  {	
	//LOAD A LOT of STUF
  	//FIRST ROUND WE NEED 0..3,128..131 LEFT SIDE (Real and Imag)
  	//256..259,384..387 (Real and Imag)
  	//REAL PART
  	v0_r  = _mm256_load_pd(x->real+offset);
  	v128_r  = _mm256_load_pd(x->real+offset+128);
  	v256_r = _mm256_load_pd(x->real+offset+256);
  	v384_r = _mm256_load_pd(x->real+offset+384);
  	//IMAG PART
	v0_i  = _mm256_load_pd(x->imag+offset);
  	v128_i  = _mm256_load_pd(x->imag+offset+128);
  	v256_i = _mm256_load_pd(x->imag+offset+256);
  	v384_i = _mm256_load_pd(x->imag+offset+384);
	  	
  	//START TWIDDLE
  	//WE ARE NOW IN 512 SO EVERYTHING BETWEEN 256 and 512 needs to be multiplied by root of unity
  	//TWIDDLE 0 - 256
	temp_imag = _mm256_mul_pd(v256_i,imag_twid_L1);	
	temp_real = _mm256_mul_pd(v256_r,imag_twid_L1);

	sub_real = _mm256_sub_pd(temp_real,temp_imag);
	temp_imag = _mm256_add_pd(temp_real,temp_imag);

	v256_r = _mm256_sub_pd(v0_r,sub_real);
	v256_i = _mm256_sub_pd(v0_i,temp_imag);

	v0_r = _mm256_add_pd(v0_r,sub_real);
	v0_i = _mm256_add_pd(v0_i,temp_imag);
	//TWIDDLE 128-384
	temp_imag = _mm256_mul_pd(v384_i,imag_twid_L1);	
	temp_real = _mm256_mul_pd(v384_r,imag_twid_L1);

	sub_real = _mm256_sub_pd(temp_real,temp_imag);
	temp_imag = _mm256_add_pd(temp_real,temp_imag);

	v384_r = _mm256_sub_pd(v128_r,sub_real);
	v384_i = _mm256_sub_pd(v128_i,temp_imag);

	v128_r = _mm256_add_pd(v128_r,sub_real);
	v128_i = _mm256_add_pd(v128_i,temp_imag);

	//NOW WE START WITH 256
	//FIRST LEFT SIDE
	//TWIDDLE 0 - 128
	temp_real = _mm256_mul_pd(v128_i,imag_twid_L2_1);
	temp_imag = _mm256_mul_pd(v128_i,real_twid_L2_1);

	temp_real = _mm256_fmsub_pd(v128_r,real_twid_L2_1,temp_real);
	temp_imag = _mm256_fmadd_pd(v128_r,imag_twid_L2_1,temp_imag);

	v128_r = _mm256_sub_pd(v0_r,temp_real);
	v128_i = _mm256_sub_pd(v0_i,temp_imag);

	v0_r = _mm256_add_pd(v0_r,temp_real);
	v0_i = _mm256_add_pd(v0_i,temp_imag);
	//NOW WE DO THE RIGHT SIDE
	//TWIDDLE 256 - 384
	temp_real = _mm256_mul_pd(v384_i,imag_twid_L2_2);
	temp_imag = _mm256_mul_pd(v384_i,real_twid_L2_2);

	temp_real = _mm256_fmsub_pd(v384_r,real_twid_L2_2,temp_real);
	temp_imag = _mm256_fmadd_pd(v384_r,imag_twid_L2_2,temp_imag);

	v384_r = _mm256_sub_pd(v256_r,temp_real);
	v384_i = _mm256_sub_pd(v256_i,temp_imag);

	v256_r = _mm256_add_pd(v256_r,temp_real);
	v256_i = _mm256_add_pd(v256_i,temp_imag);
	
	//STORE ALL RESULTS
	_mm256_store_pd(x->real+offset,v0_r);
	_mm256_store_pd(x->real+offset+128,v128_r);
	_mm256_store_pd(x->real+offset+256,v256_r);
	_mm256_store_pd(x->real+offset+384,v384_r);

	_mm256_store_pd(x->imag+offset,v0_i);
	_mm256_store_pd(x->imag+offset+128,v128_i);
	_mm256_store_pd(x->imag+offset+256,v256_i);
	_mm256_store_pd(x->imag+offset+384,v384_i);
 }

  //Merge 128-64
  for(lo = 0;lo < CPLXDIM;lo +=128)
  {	
	//init twiddle for 128
	real_twid_L1 = _mm256_set1_pd(wortel[0][6][lo/128]);
	imag_twid_L1 = _mm256_set1_pd(wortel[1][6][lo/128]);
	//init twiddle for 64
	real_twid_L2_1 = _mm256_set1_pd(wortel[0][5][lo/64]);
	imag_twid_L2_1 = _mm256_set1_pd(wortel[1][5][lo/64]);
	real_twid_L2_2 = _mm256_set1_pd(wortel[0][5][(lo+64)/64]);
	imag_twid_L2_2 = _mm256_set1_pd(wortel[1][5][(lo+64)/64]);

  	for(int offset = 0; offset < 32; offset+=4)
  	{	
  		//LOAD A LOT of STUF
	  	//FIRST ROUND WE NEED 0..3,32..35 LEFT SIDE (Real and Imag)
	  	//64..67,96..99 (Real and Imag)
	  	//REAL PART
	  	v0_r   = _mm256_load_pd(x->real+lo+offset);
	  	v128_r = _mm256_load_pd(x->real+lo+offset+32);
	  	v256_r = _mm256_load_pd(x->real+lo+offset+64);
	  	v384_r = _mm256_load_pd(x->real+lo+offset+96);
	  	//IMAG PART
		v0_i  = _mm256_load_pd(x->imag+lo+offset);
	  	v128_i  = _mm256_load_pd(x->imag+lo+offset+32);
	  	v256_i = _mm256_load_pd(x->imag+lo+offset+64);
	  	v384_i = _mm256_load_pd(x->imag+lo+offset+96);
	  	
	  	//TWIDDLE 0 - 64
		temp_real = _mm256_mul_pd(v256_i,imag_twid_L1);
		temp_imag = _mm256_mul_pd(v256_i,real_twid_L1);

		temp_real = _mm256_fmsub_pd(v256_r,real_twid_L1,temp_real);
		temp_imag = _mm256_fmadd_pd(v256_r,imag_twid_L1,temp_imag);

		v256_r = _mm256_sub_pd(v0_r,temp_real);
		v256_i = _mm256_sub_pd(v0_i,temp_imag);

		v0_r = _mm256_add_pd(v0_r,temp_real);
		v0_i = _mm256_add_pd(v0_i,temp_imag);
	  	//TWIDDLE 32 - 96
		temp_real = _mm256_mul_pd(v384_i,imag_twid_L1);
		temp_imag = _mm256_mul_pd(v384_i,real_twid_L1);

		temp_real = _mm256_fmsub_pd(v384_r,real_twid_L1,temp_real);
		temp_imag = _mm256_fmadd_pd(v384_r,imag_twid_L1,temp_imag);

		v384_r = _mm256_sub_pd(v128_r,temp_real);
		v384_i = _mm256_sub_pd(v128_i,temp_imag);

		v128_r = _mm256_add_pd(v128_r,temp_real);
		v128_i = _mm256_add_pd(v128_i,temp_imag);

		//NOW WE START WITH 64
		//FIRST LEFT SIDE
		//TWIDDLE 0 - 32
		temp_real = _mm256_mul_pd(v128_i,imag_twid_L2_1);
		temp_imag = _mm256_mul_pd(v128_i,real_twid_L2_1);

		temp_real = _mm256_fmsub_pd(v128_r,real_twid_L2_1,temp_real);
		temp_imag = _mm256_fmadd_pd(v128_r,imag_twid_L2_1,temp_imag);

		v128_r = _mm256_sub_pd(v0_r,temp_real);
		v128_i = _mm256_sub_pd(v0_i,temp_imag);

		v0_r = _mm256_add_pd(v0_r,temp_real);
		v0_i = _mm256_add_pd(v0_i,temp_imag);
		//NOW WE DO THE RIGHT SIDE
		//TWIDDLE 64 - 96
		temp_real = _mm256_mul_pd(v384_i,imag_twid_L2_2);
		temp_imag = _mm256_mul_pd(v384_i,real_twid_L2_2);

		temp_real = _mm256_fmsub_pd(v384_r,real_twid_L2_2,temp_real);
		temp_imag = _mm256_fmadd_pd(v384_r,imag_twid_L2_2,temp_imag);

		v384_r = _mm256_sub_pd(v256_r,temp_real);
		v384_i = _mm256_sub_pd(v256_i,temp_imag);

		v256_r = _mm256_add_pd(v256_r,temp_real);
		v256_i = _mm256_add_pd(v256_i,temp_imag);

		//STORE ALL RESULTS
		_mm256_store_pd(x->real+lo+offset,v0_r);
		_mm256_store_pd(x->real+lo+offset+32,v128_r);
		_mm256_store_pd(x->real+lo+offset+64,v256_r);
		_mm256_store_pd(x->real+lo+offset+96,v384_r);

		_mm256_store_pd(x->imag+lo+offset,v0_i);
		_mm256_store_pd(x->imag+lo+offset+32,v128_i);
		_mm256_store_pd(x->imag+lo+offset+64,v256_i);
		_mm256_store_pd(x->imag+lo+offset+96,v384_i);
  	}
  }

    //Merge 32-16
  for(lo = 0;lo < CPLXDIM;lo +=32)
  {	
	//init twiddle for 32
	real_twid_L1 = _mm256_set1_pd(wortel[0][4][lo/32]);
	imag_twid_L1 = _mm256_set1_pd(wortel[1][4][lo/32]);
	//init twiddle for 16
	real_twid_L2_1 = _mm256_set1_pd(wortel[0][3][lo/16]);
	imag_twid_L2_1 = _mm256_set1_pd(wortel[1][3][lo/16]);
	real_twid_L2_2 = _mm256_set1_pd(wortel[0][3][(lo+16)/16]);
	imag_twid_L2_2 = _mm256_set1_pd(wortel[1][3][(lo+16)/16]);

  	for(int offset = 0; offset < 8; offset+=4)
  	{	
  		//LOAD A LOT of STUF
	  	//FIRST ROUND WE NEED 0..3,32..35 LEFT SIDE (Real and Imag)
	  	//64..67,96..99 (Real and Imag)
	  	//REAL PART
	  	v0_r   = _mm256_load_pd(x->real+lo+offset);
	  	v128_r = _mm256_load_pd(x->real+lo+offset+8);
	  	v256_r = _mm256_load_pd(x->real+lo+offset+16);
	  	v384_r = _mm256_load_pd(x->real+lo+offset+24);
	  	//IMAG PART
		v0_i  = _mm256_load_pd(x->imag+lo+offset);
	  	v128_i  = _mm256_load_pd(x->imag+lo+offset+8);
	  	v256_i = _mm256_load_pd(x->imag+lo+offset+16);
	  	v384_i = _mm256_load_pd(x->imag+lo+offset+24);
	  	
		//TWIDDLE 0 - 16
		temp_real = _mm256_mul_pd(v256_i,imag_twid_L1);
		temp_imag = _mm256_mul_pd(v256_i,real_twid_L1);

		temp_real = _mm256_fmsub_pd(v256_r,real_twid_L1,temp_real);
		temp_imag = _mm256_fmadd_pd(v256_r,imag_twid_L1,temp_imag);

		v256_r = _mm256_sub_pd(v0_r,temp_real);
		v256_i = _mm256_sub_pd(v0_i,temp_imag);

		v0_r = _mm256_add_pd(v0_r,temp_real);
		v0_i = _mm256_add_pd(v0_i,temp_imag);
	  	//TWIDDLE 8 - 24
		temp_real = _mm256_mul_pd(v384_i,imag_twid_L1);
		temp_imag = _mm256_mul_pd(v384_i,real_twid_L1);

		temp_real = _mm256_fmsub_pd(v384_r,real_twid_L1,temp_real);
		temp_imag = _mm256_fmadd_pd(v384_r,imag_twid_L1,temp_imag);

		v384_r = _mm256_sub_pd(v128_r,temp_real);
		v384_i = _mm256_sub_pd(v128_i,temp_imag);

		v128_r = _mm256_add_pd(v128_r,temp_real);
		v128_i = _mm256_add_pd(v128_i,temp_imag);

		//NOW WE START WITH 16
		//FIRST LEFT SIDE
		//TWIDDLE 0 - 8
		temp_real = _mm256_mul_pd(v128_i,imag_twid_L2_1);
		temp_imag = _mm256_mul_pd(v128_i,real_twid_L2_1);

		temp_real = _mm256_fmsub_pd(v128_r,real_twid_L2_1,temp_real);
		temp_imag = _mm256_fmadd_pd(v128_r,imag_twid_L2_1,temp_imag);

		v128_r = _mm256_sub_pd(v0_r,temp_real);
		v128_i = _mm256_sub_pd(v0_i,temp_imag);

		v0_r = _mm256_add_pd(v0_r,temp_real);
		v0_i = _mm256_add_pd(v0_i,temp_imag);
		//NOW WE DO THE RIGHT SIDE
		//TWIDDLE 16 - 24
		temp_real = _mm256_mul_pd(v384_i,imag_twid_L2_2);
		temp_imag = _mm256_mul_pd(v384_i,real_twid_L2_2);

		temp_real = _mm256_fmsub_pd(v384_r,real_twid_L2_2,temp_real);
		temp_imag = _mm256_fmadd_pd(v384_r,imag_twid_L2_2,temp_imag);

		v384_r = _mm256_sub_pd(v256_r,temp_real);
		v384_i = _mm256_sub_pd(v256_i,temp_imag);

		v256_r = _mm256_add_pd(v256_r,temp_real);
		v256_i = _mm256_add_pd(v256_i,temp_imag);

		//STORE ALL RESULTS
		_mm256_store_pd(x->real+lo+offset,v0_r);
		_mm256_store_pd(x->real+lo+offset+8,v128_r);
		_mm256_store_pd(x->real+lo+offset+16,v256_r);
		_mm256_store_pd(x->real+lo+offset+24,v384_r);

		_mm256_store_pd(x->imag+lo+offset,v0_i);
		_mm256_store_pd(x->imag+lo+offset+8,v128_i);
		_mm256_store_pd(x->imag+lo+offset+16,v256_i);
		_mm256_store_pd(x->imag+lo+offset+24,v384_i);
  	}
  }

  //Do layer 8-4-2
  for(lo=0;lo < CPLXDIM;lo +=8)
  {	
    
    real_twid_L1 = _mm256_set1_pd(wortel[0][2][lo/8]);
    imag_twid_L1 = _mm256_set1_pd(wortel[1][2][lo/8]);
	//(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
	v0_r = _mm256_load_pd(x->real+lo);
	v0_i = _mm256_load_pd(x->imag+lo);
	v256_r = _mm256_load_pd(x->real+lo+4);
	v256_i = _mm256_load_pd(x->imag+lo+4);
	//TEMP_real = bd
	temp_real = _mm256_mul_pd(v256_i,imag_twid_L1);
	//TEMP_imag = bc
	temp_imag = _mm256_mul_pd(v256_i,real_twid_L1);

	//TEMP_real = ac - bd
	temp_real = _mm256_fmsub_pd(v256_r,real_twid_L1,temp_real);
	//TEMP_imag = ad + bc
	temp_imag = _mm256_fmadd_pd(v256_r,imag_twid_L1,temp_imag);

	v256_r = _mm256_sub_pd(v0_r,temp_real);
	v256_i = _mm256_sub_pd(v0_i,temp_imag);

	v0_r = _mm256_add_pd(v0_r,temp_real);
	v0_i = _mm256_add_pd(v0_i,temp_imag);

	//START LAYER 4
	//We know that v0_r = |a|b|c|d| v0_i = |ai|bi|ci|di| v256_r = |e|f|g|h| v256_i = |ei|fi|gi|hi| 
	//We need to twidle (c+ci),(d+di) and (g+gi),(h+hi)

	sub_real = _mm256_permute2f128_pd(v0_r,v256_r,0x31);
	sub_imag = _mm256_permute2f128_pd(v0_i,v256_i,0x31);
	real_twid_L1 = _mm256_setr_pd(wortel[0][1][lo/4],wortel[0][1][lo/4],wortel[0][1][(lo+4)/4],wortel[0][1][(lo+4)/4]);
    imag_twid_L1 = _mm256_setr_pd(wortel[1][1][lo/4],wortel[1][1][lo/4],wortel[1][1][(lo+4)/4],wortel[1][1][(lo+4)/4]);

    temp_real = _mm256_mul_pd(sub_imag,imag_twid_L1);
    temp_imag = _mm256_mul_pd(sub_imag,real_twid_L1);
	//TEMP_real = ac - bd
	temp_real = _mm256_fmsub_pd(sub_real,real_twid_L1,temp_real);
	//TEMP_imag = ad + bc
	temp_imag = _mm256_fmadd_pd(sub_real,imag_twid_L1,temp_imag);
	
	//REAL PART
	//get abef
	sub_real = _mm256_permute2f128_pd(v0_r,v256_r,0x20);
	//abef  
	//cdgh-
	sub_imag = _mm256_sub_pd(sub_real,temp_real);
	//abef
	//cdgh+
	sub_real = _mm256_add_pd(sub_real,temp_real);

	// //NEEDED TO COMPLETE LAYER 4
	// v0_r = _mm256_permute2f128_pd(sub_imag,sub_real,0x02);
	// v256_r = _mm256_permute2f128_pd(sub_imag,sub_real,0x13);
	//PREPARE REALS FOR LAYER 2
	//STORE ALL FACTORS THAT NEED TO BE MULT WITH ROOTS OF UNITY IN v0_r
	v0_r = _mm256_unpackhi_pd(sub_real,sub_imag);
	v256_r = _mm256_unpacklo_pd(sub_real,sub_imag);


	//IMAG PART
	//get ai bi ei fi
	sub_real = _mm256_permute2f128_pd(v0_i,v256_i,0x20);
	//abef  
	//cdgh-
	sub_imag = _mm256_sub_pd(sub_real,temp_imag);
	//abef
	//cdgh+
	sub_real = _mm256_add_pd(sub_real,temp_imag);

	//NEEDED TO COMPLETE LAYER 4
	// v0_i = _mm256_permute2f128_pd(sub_imag,sub_real,0x02);
	// v256_i = _mm256_permute2f128_pd(sub_imag,sub_real,0x13);
	//PREPARE IMAGS FOR LAYER 2
	//STORE ALL FACTORS THAT NEED TO BE MULT WITH ROOTS OF UNITY IN v0_i
	//STORE ALL NON TWIDLE IMAGS IN v256_i 
	v0_i = _mm256_unpackhi_pd(sub_real,sub_imag);
	v256_i = _mm256_unpacklo_pd(sub_real,sub_imag);

	// // //START LAYER 2!!
	real_twid_L1 = _mm256_setr_pd(wortel[0][0][lo/2],wortel[0][0][(lo+2)/2],wortel[0][0][(lo+4)/2],wortel[0][0][(lo+6)/2]);
    imag_twid_L1 = _mm256_setr_pd(wortel[1][0][lo/2],wortel[1][0][(lo+2)/2],wortel[1][0][(lo+4)/2],wortel[1][0][(lo+6)/2]);

    temp_real = _mm256_mul_pd(v0_i,imag_twid_L1);
    temp_imag = _mm256_mul_pd(v0_i,real_twid_L1);

	//TEMP_real = ac - bd
	temp_real = _mm256_fmsub_pd(v0_r,real_twid_L1,temp_real);
	//TEMP_imag = ad + bc
	temp_imag = _mm256_fmadd_pd(v0_r,imag_twid_L1,temp_imag);

	sub_real = _mm256_sub_pd(v256_r,temp_real);
	sub_imag = _mm256_add_pd(v256_r,temp_real); 

	temp_real = _mm256_unpacklo_pd(sub_imag,sub_real);
	sub_real  = _mm256_unpackhi_pd(sub_imag,sub_real);

	v0_r = _mm256_permute2f128_pd(temp_real,sub_real,0x20);
	v256_r = _mm256_permute2f128_pd(temp_real,sub_real,0x31);

	sub_real = _mm256_sub_pd(v256_i,temp_imag);
	sub_imag = _mm256_add_pd(v256_i,temp_imag);

	temp_real = _mm256_unpacklo_pd(sub_imag,sub_real);
	sub_real  = _mm256_unpackhi_pd(sub_imag,sub_real);

	v0_i = _mm256_permute2f128_pd(temp_real,sub_real,0x20);
	v256_i = _mm256_permute2f128_pd(temp_real,sub_real,0x31);

	_mm256_store_pd(x->real+lo,v0_r);
	_mm256_store_pd(x->imag+lo,v0_i);
	_mm256_store_pd(x->real+lo+4,v256_r);
	_mm256_store_pd(x->imag+lo+4,v256_i);	
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
  iterative_phi(x);
}

void phi_backward(cplx_ptr *x, ring_t *ring)
{
  inverse_phi(x);

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    ring->v[i] = x->real[i];
    ring->v[j] = x->imag[i];
    ++j;
  }
}